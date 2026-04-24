#!/usr/bin/env python3
"""
CASP17 Spider - Recursively crawls CASP17 pages from predictioncenter.org

CASP17 hasn't started yet, so this captures all announcement/info pages
and converts them to readable markdown. Strict CASP17 scope — URLs
outside /casp17/ or /download_area/CASP17/ are ignored.

Usage:
    python3 spider17.py --dry-run
    python3 spider17.py
    python3 spider17.py --resume
    python3 spider17.py --max-pages 1000
"""

import argparse
import json
import re
import time
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from urllib.parse import urljoin, urlparse, urldefrag, parse_qsl, urlencode

import requests
from bs4 import BeautifulSoup
from markdownify import markdownify as md

BASE_URL = "https://predictioncenter.org"
HOST = "predictioncenter.org"
OUTPUT_DIR = Path(__file__).parent / "casp17"
STATE_FILE = Path(__file__).parent / ".spider17_state.json"
DEFAULT_DELAY = 1.0
MAX_DOWNLOAD_BYTES = 20 * 1024 * 1024

SEED_URLS = [
    f"{BASE_URL}/casp17/",
    f"{BASE_URL}/casp17/index.html",
    f"{BASE_URL}/casp17/index.cgi",
    f"{BASE_URL}/download_area/CASP17/",
]

# Binary/media we don't want to parse as HTML but still worth keeping
DOWNLOADABLE_EXT = {
    ".txt", ".csv", ".tsv", ".seq", ".fasta", ".pdb", ".cif",
    ".json", ".xml", ".yaml", ".yml", ".log",
}

# Binary we skip entirely (too large or not useful as text)
SKIP_EXT = {
    ".jpg", ".jpeg", ".png", ".gif", ".svg", ".ico",
    ".mp3", ".mp4", ".mov", ".avi",
    ".zip", ".tar", ".gz", ".tgz", ".bz2", ".xz", ".7z",
    ".pdf", ".ps", ".doc", ".docx",
}

# Other CASP rounds to explicitly avoid
OTHER_CASP_RE = re.compile(r"/casp(?:1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16)(?:[/_-]|$)", re.I)
CASP17_RE = re.compile(r"casp17", re.I)

# Query params to strip for dedup (sort controls on targetlist.cgi etc.)
STRIP_PARAMS = {"order", "field"}


def normalize(url: str) -> str:
    url = urldefrag(url).url.rstrip("#")
    parsed = urlparse(url)
    if parsed.query:
        params = sorted(
            (k, v) for k, v in parse_qsl(parsed.query, keep_blank_values=True)
            if k not in STRIP_PARAMS
        )
        parsed = parsed._replace(query=urlencode(params))
        url = parsed.geturl()
    return url


def in_scope(url: str) -> bool:
    parsed = urlparse(url)
    if parsed.scheme not in ("http", "https"):
        return False
    if parsed.netloc and parsed.netloc != HOST:
        return False
    path_and_query = (parsed.path or "") + ("?" + parsed.query if parsed.query else "")
    if not CASP17_RE.search(path_and_query):
        return False
    if OTHER_CASP_RE.search(parsed.path or ""):
        return False
    return True


def url_to_path(url: str) -> Path:
    parsed = urlparse(url)
    path = parsed.path or "/"
    if path.endswith("/"):
        path += "index.html"

    # Map known prefixes to clean output dirs
    if re.match(r"^/casp17/?", path, re.I):
        rel = re.sub(r"^/casp17/?", "", path, flags=re.I)
    elif re.match(r"^/download_area/CASP17/?", path, re.I):
        rel = "downloads/" + re.sub(r"^/download_area/CASP17/?", "", path, flags=re.I)
    else:
        rel = path.lstrip("/")

    if not rel:
        rel = "index.html"

    # Swap HTML/CGI extensions for .md
    rel, suffix = re.subn(r"\.(html?|cgi)$", "", rel, flags=re.I)
    ext_swapped = suffix > 0

    if parsed.query:
        safe_q = re.sub(r"[^a-zA-Z0-9._=-]", "_", parsed.query)
        rel = f"{rel}__{safe_q}"

    if ext_swapped or "." not in Path(rel).name:
        rel = rel + ".md"

    return OUTPUT_DIR / rel


@dataclass
class SpiderState:
    crawled_urls: set = field(default_factory=set)
    downloaded_files: set = field(default_factory=set)
    failed_urls: dict = field(default_factory=dict)
    queue: list = field(default_factory=list)

    def save(self, path: Path):
        data = {
            "crawled_urls": sorted(self.crawled_urls),
            "downloaded_files": sorted(self.downloaded_files),
            "failed_urls": self.failed_urls,
            "queue": self.queue,
        }
        path.write_text(json.dumps(data, indent=2))

    @classmethod
    def load(cls, path: Path) -> "SpiderState":
        if not path.exists():
            return cls()
        data = json.loads(path.read_text())
        s = cls()
        s.crawled_urls = set(data.get("crawled_urls", []))
        s.downloaded_files = set(data.get("downloaded_files", []))
        s.failed_urls = data.get("failed_urls", {})
        s.queue = data.get("queue", [])
        return s


class Casp17Spider:
    def __init__(self, dry_run=False, delay=DEFAULT_DELAY, resume=False, max_pages=500):
        self.dry_run = dry_run
        self.delay = delay
        self.max_pages = max_pages
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "CASPSpider/1.0 (Academic research; sness@sness.net)"
        })
        self.state = SpiderState.load(STATE_FILE) if resume else SpiderState()

    def log(self, msg, level="INFO"):
        prefix = "[DRY-RUN] " if self.dry_run else ""
        print(f"{prefix}[{level}] {msg}")

    def fetch(self, url: str, stream: bool = False) -> Optional[requests.Response]:
        time.sleep(self.delay)
        try:
            resp = self.session.get(url, timeout=30, stream=stream)
            resp.raise_for_status()
            return resp
        except requests.RequestException as e:
            self.log(f"Failed: {url}: {e}", "ERROR")
            self.state.failed_urls[url] = str(e)
            return None

    def extract_links(self, html: str, base_url: str) -> list:
        soup = BeautifulSoup(html, "html.parser")
        out = []
        for a in soup.find_all("a", href=True):
            href = a["href"].strip()
            if not href or href.startswith(("mailto:", "javascript:", "#")):
                continue
            absolute = normalize(urljoin(base_url, href))
            out.append(absolute)
        return out

    @staticmethod
    def _is_data_table(table) -> bool:
        rows = table.find_all("tr")
        if len(rows) < 2:
            return False
        for row in rows:
            for cell in row.find_all(["td", "th"]):
                if cell.find(["table", "div", "h1", "h2", "h3", "form", "ul", "ol", "p"]):
                    return False
        if table.find("th") and len(rows) >= 2:
            return True
        return len(rows) >= 3

    @staticmethod
    def _flatten_table(table):
        for tag in [table] + list(table.find_all(["tr", "td", "th", "tbody", "thead", "tfoot"])):
            tag.name = "div"
            tag.attrs = {}

    def _process_tables(self, soup):
        data_tags = ("_dtable", "_dtr", "_dtd", "_dth")
        while True:
            ready = [t for t in soup.find_all("table") if not t.find("table")]
            if not ready:
                break
            for table in ready:
                if self._is_data_table(table):
                    table.name = data_tags[0]
                    for tr in table.find_all("tr"):
                        tr.name = data_tags[1]
                    for td in table.find_all("td"):
                        td.name = data_tags[2]
                    for th in table.find_all("th"):
                        th.name = data_tags[3]
                else:
                    self._flatten_table(table)
        for original, dtag in zip(("table", "tr", "td", "th"), data_tags):
            for t in soup.find_all(dtag):
                t.name = original

    def html_to_markdown(self, html: str, url: str) -> str:
        soup = BeautifulSoup(html, "html.parser")
        for el in soup(["script", "style", "noscript"]):
            el.decompose()
        for el in soup.find_all("div", class_="dtree"):
            el.decompose()
        for el in soup.find_all("div", id="treemenu"):
            el.decompose()

        self._process_tables(soup)

        content = (
            soup.find("div", id="content")
            or soup.find("div", class_="content")
            or soup.find("main")
            or soup.body
            or soup
        )

        title_el = soup.find("title")
        title = re.sub(r"\s+", " ", title_el.get_text().strip()) if title_el else url

        body = md(str(content), heading_style="ATX").strip()
        body = re.sub(r"\n{3,}", "\n\n", body)
        body = re.sub(r"[ \t]+\n", "\n", body)
        body = re.sub(r"(?m)^\s*\|\s*\|\s*$", "", body)

        header = (
            f"<!-- Source: {url} -->\n"
            f"<!-- Crawled: {time.strftime('%Y-%m-%d %H:%M:%S')} -->\n\n"
            f"# {title}\n\n"
        )
        return header + body + "\n"

    def save_page(self, url: str, html: str) -> Path:
        output = url_to_path(url)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(self.html_to_markdown(html, url))
        return output

    def download_binary(self, url: str) -> Optional[Path]:
        resp = self.fetch(url, stream=True)
        if resp is None:
            return None

        output = url_to_path(url)
        ext = Path(urlparse(url).path).suffix
        if ext and output.suffix == ".md":
            output = output.with_suffix(ext)

        output.parent.mkdir(parents=True, exist_ok=True)
        size = 0
        try:
            with open(output, "wb") as f:
                for chunk in resp.iter_content(chunk_size=16384):
                    size += len(chunk)
                    if size > MAX_DOWNLOAD_BYTES:
                        f.close()
                        output.unlink(missing_ok=True)
                        self.log(f"Aborted (>{MAX_DOWNLOAD_BYTES} bytes): {url}", "WARN")
                        return None
                    f.write(chunk)
        except requests.RequestException as e:
            self.log(f"Download failed: {url}: {e}", "ERROR")
            self.state.failed_urls[url] = str(e)
            return None

        self.state.downloaded_files.add(str(output))
        self.log(f"Downloaded: {output} ({size} bytes)")
        return output

    def run(self):
        self.log(f"Starting CASP17 spider (dry_run={self.dry_run})")

        queue = deque(self.state.queue) if self.state.queue else deque(SEED_URLS)
        seen_in_queue = set(queue)
        visited = 0

        try:
            while queue and visited < self.max_pages:
                url = normalize(queue.popleft())
                if url in self.state.crawled_urls:
                    continue
                if not in_scope(url):
                    continue

                ext = Path(urlparse(url).path).suffix.lower()
                if ext in SKIP_EXT:
                    self.state.crawled_urls.add(url)
                    continue

                if self.dry_run:
                    self.log(f"Would fetch: {url}")
                    self.state.crawled_urls.add(url)
                    visited += 1
                    continue

                looks_like_page = ext in ("", ".html", ".htm", ".cgi")
                if not looks_like_page and ext in DOWNLOADABLE_EXT:
                    self.download_binary(url)
                    self.state.crawled_urls.add(url)
                    visited += 1
                    continue

                resp = self.fetch(url)
                if resp is None:
                    self.state.crawled_urls.add(url)
                    continue

                self.state.crawled_urls.add(url)
                visited += 1

                content_type = resp.headers.get("content-type", "").lower()
                if "html" in content_type or looks_like_page:
                    html = resp.text
                    output = self.save_page(url, html)
                    self.log(f"Saved: {output}")

                    for link in self.extract_links(html, url):
                        if link in self.state.crawled_urls or link in seen_in_queue:
                            continue
                        if not in_scope(link):
                            continue
                        queue.append(link)
                        seen_in_queue.add(link)
                elif "text" in content_type:
                    output = url_to_path(url)
                    if output.suffix == ".md" and Path(urlparse(url).path).suffix:
                        output = output.with_suffix(Path(urlparse(url).path).suffix)
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_text(resp.text)
                    self.log(f"Saved text: {output}")
                else:
                    self.log(f"Skipping unknown content-type {content_type}: {url}")

        except KeyboardInterrupt:
            self.log("Interrupted by user", "WARN")

        self.state.queue = list(queue)
        if not self.dry_run:
            self.state.save(STATE_FILE)
            self.log(f"State saved to {STATE_FILE}")

        self.log("=== Spider Complete ===")
        self.log(f"Visited: {visited}  Crawled total: {len(self.state.crawled_urls)}")
        self.log(f"Downloaded: {len(self.state.downloaded_files)} files")
        self.log(f"Queue remaining: {len(self.state.queue)}")
        if self.state.failed_urls:
            self.log(f"Failed: {len(self.state.failed_urls)} URLs", "WARN")


def main():
    parser = argparse.ArgumentParser(description="Recursively spider CASP17 pages")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--delay", type=float, default=DEFAULT_DELAY)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max-pages", type=int, default=500)
    args = parser.parse_args()

    spider = Casp17Spider(
        dry_run=args.dry_run,
        delay=args.delay,
        resume=args.resume,
        max_pages=args.max_pages,
    )
    spider.run()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
CASP Spider - Crawls CASP15 and CASP16 ligand prediction data from predictioncenter.org

Usage:
    python spider.py --dry-run           # Show what would be crawled
    python spider.py                     # Run full crawl
    python spider.py --section ligand    # Only crawl ligand sections
    python spider.py --resume            # Resume from saved state
"""

import argparse
import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup

# Configuration
BASE_URL = "https://predictioncenter.org"
OUTPUT_DIR = Path(__file__).parent
STATE_FILE = OUTPUT_DIR / ".spider_state.json"
MAX_FILE_SIZE = 1 * 1024 * 1024  # 1MB limit for downloads
DEFAULT_DELAY = 1.0  # seconds between requests

# CASP16 ligand targets
CASP16_LIGAND_TARGETS = ["D1273", "R1261v1", "R1262v1", "R1263v1", "R1264v1", "R1288", "T1214"]

# CASP15 ligand targets
CASP15_LIGAND_TARGETS = [
    "H1114", "H1135", "H1171v1", "H1171v2",
    "H1172v1", "H1172v2", "H1172v3", "H1172v4",
    "R1117v2", "T1118v1", "T1124", "T1127v2",
    "T1146", "T1152",
    "T1158v1", "T1158v2", "T1158v3", "T1158v4",
    "T1170", "T1181", "T1186", "T1187", "T1188"
]

# Files to download from pharma_ligands directory
PHARMA_LIGAND_FILES = [
    "L1000_chymase_README.txt",
    "L1000_exper_affinity.csv",
    "L1000_exper_affinity_obsolete.csv",
    "L1000.SMILES.tar.gz",
    "L2000_cathepsin_README.txt",
    "L2000.SMILES.tar.gz",
    "L3000_autotaxin_README.txt",
    "L3000_exper_affinity.csv",
    "L3000.SMILES.tar.gz",
    "L4000_mpro_README.txt",
    "L4000.SMILES.tar.gz",
    "L4020_updated.tsv",
]

# Sequence files
CASP16_SEQUENCE_FILES = [
    "casp16.H1.seq.txt",
    "casp16.RDM1.seq.txt",
    "casp16.T1.seq.txt",
]

# Ligand results/scoring files
LIGAND_RESULTS_FILES = [
    # CASP16 scoring results
    ("CASP16/results/ligands/RESULTS.pose.SUMMARY.csv", "casp16/results/ligand_scores.csv", 15 * 1024 * 1024),
    ("CASP16/results/ligands/RESULTS.pose.L5001v1.csv", "casp16/results/ligand_scores_L5001v1.csv", 1 * 1024 * 1024),
    # CASP15 scoring results
    ("CASP15/results/tables/ligand.csv", "casp15/results/ligand_scores.csv", 10 * 1024 * 1024),
]

# CASP16 ligand prediction archives (non-pharma, smaller files)
CASP16_LIGAND_PREDICTIONS = [
    ("D1273.tar.gz", 2 * 1024 * 1024),      # 1.2M
    ("R1261.tar.gz", 6 * 1024 * 1024),      # 5.3M
    ("R1261v1.tgz", 4 * 1024 * 1024),       # 3.2M
    ("R1262.tar.gz", 6 * 1024 * 1024),      # 5.1M
    ("R1262v1.tgz", 4 * 1024 * 1024),       # 3.2M
    ("R1263.tar.gz", 4 * 1024 * 1024),      # 3.5M
    ("R1263v1.tgz", 3 * 1024 * 1024),       # 2.4M
    ("R1264.tar.gz", 4 * 1024 * 1024),      # 3.6M
    ("R1264v1.tgz", 3 * 1024 * 1024),       # 2.4M
    ("R1288.tar.gz", 3 * 1024 * 1024),      # 2.7M
    ("T1214.tar.gz", 12 * 1024 * 1024),     # 10M
    ("L5001v1.tgz", 12 * 1024 * 1024),      # 9.7M
]

# CASP16 pharma affinity files (small)
CASP16_PHARMA_AFFINITY = [
    ("L1000.AFFNTY", 100 * 1024),           # 63K
    ("L3000.AFFNTY", 500 * 1024),           # 377K
    ("STAGE2.AFFNTY.tgz", 50 * 1024),       # 27K
]

# CASP15 ligand prediction archives
CASP15_LIGAND_PREDICTIONS = [
    ("R1117.tar.gz", 1 * 1024 * 1024),      # 646K
    ("R1117v2.tar.gz", 1 * 1024 * 1024),    # 674K
    ("T1152.tar.gz", 3 * 1024 * 1024),      # 2.1M
    ("T1127.tar.gz", 4 * 1024 * 1024),      # 3.4M
    ("T1186.tar.gz", 4 * 1024 * 1024),      # 3.6M
    ("R1126.tar.gz", 6 * 1024 * 1024),      # 4.9M
    ("T1105v1.tar.gz", 6 * 1024 * 1024),    # 5.2M
    ("T1146.tar.gz", 6 * 1024 * 1024),      # 5.3M
    ("T1187.tar.gz", 6 * 1024 * 1024),      # 5.4M
    ("T1127v2.tar.gz", 7 * 1024 * 1024),    # 5.8M
    ("R1136.tar.gz", 8 * 1024 * 1024),      # 7.0M
    ("T1188.tar.gz", 10 * 1024 * 1024),     # 8.8M
    ("T1118v1.tar.gz", 12 * 1024 * 1024),   # 10M
    ("T1118.tar.gz", 12 * 1024 * 1024),     # 11M
    ("T1124.tar.gz", 14 * 1024 * 1024),     # 12M
    ("T1158v1.tar.gz", 22 * 1024 * 1024),   # 19M
    ("T1158v2.tar.gz", 22 * 1024 * 1024),   # 19M
    ("T1158v3.tar.gz", 22 * 1024 * 1024),   # 19M
    ("T1158v4.tar.gz", 22 * 1024 * 1024),   # 19M
    ("H1135.tar.gz", 25 * 1024 * 1024),     # 20M
    ("T1170.tar.gz", 25 * 1024 * 1024),     # 21M
    ("H1171.tar.gz", 30 * 1024 * 1024),     # 26M
    ("H1172.tar.gz", 30 * 1024 * 1024),     # 26M
    ("T1181.tar.gz", 30 * 1024 * 1024),     # 27M
    ("H1114.tar.gz", 80 * 1024 * 1024),     # 76M
]


@dataclass
class SpiderState:
    """Tracks crawl progress for resumption"""
    crawled_urls: set = field(default_factory=set)
    downloaded_files: set = field(default_factory=set)
    failed_urls: dict = field(default_factory=dict)

    def save(self, path: Path):
        data = {
            "crawled_urls": list(self.crawled_urls),
            "downloaded_files": list(self.downloaded_files),
            "failed_urls": self.failed_urls,
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, path: Path) -> "SpiderState":
        if not path.exists():
            return cls()
        with open(path) as f:
            data = json.load(f)
        state = cls()
        state.crawled_urls = set(data.get("crawled_urls", []))
        state.downloaded_files = set(data.get("downloaded_files", []))
        state.failed_urls = data.get("failed_urls", {})
        return state


class CASPSpider:
    def __init__(self, dry_run: bool = False, delay: float = DEFAULT_DELAY, resume: bool = False):
        self.dry_run = dry_run
        self.delay = delay
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "CASPSpider/1.0 (Academic research; contact@example.com)"
        })
        self.state = SpiderState.load(STATE_FILE) if resume else SpiderState()

    def log(self, msg: str, level: str = "INFO"):
        prefix = "[DRY-RUN] " if self.dry_run else ""
        print(f"{prefix}[{level}] {msg}")

    def fetch(self, url: str) -> Optional[requests.Response]:
        """Fetch a URL with rate limiting"""
        if url in self.state.crawled_urls:
            self.log(f"Skipping (already crawled): {url}")
            return None

        if self.dry_run:
            self.log(f"Would fetch: {url}")
            return None

        self.log(f"Fetching: {url}")
        time.sleep(self.delay)

        try:
            resp = self.session.get(url, timeout=30)
            resp.raise_for_status()
            self.state.crawled_urls.add(url)
            return resp
        except requests.RequestException as e:
            self.log(f"Failed to fetch {url}: {e}", "ERROR")
            self.state.failed_urls[url] = str(e)
            return None

    def download_file(self, url: str, output_path: Path, max_size: int = MAX_FILE_SIZE) -> bool:
        """Download a file with size limit"""
        if str(output_path) in self.state.downloaded_files:
            self.log(f"Skipping (already downloaded): {output_path}")
            return True

        if self.dry_run:
            self.log(f"Would download: {url} -> {output_path}")
            return True

        self.log(f"Downloading: {url}")
        time.sleep(self.delay)

        try:
            # First check content-length
            resp = self.session.head(url, timeout=10)
            content_length = int(resp.headers.get("content-length", 0))
            if content_length > max_size:
                self.log(f"Skipping (too large: {content_length} bytes): {url}", "WARN")
                return False

            # Download
            resp = self.session.get(url, timeout=60, stream=True)
            resp.raise_for_status()

            # Check actual size while downloading
            output_path.parent.mkdir(parents=True, exist_ok=True)
            size = 0
            with open(output_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    size += len(chunk)
                    if size > max_size:
                        self.log(f"Aborting (exceeded size limit): {url}", "WARN")
                        f.close()
                        output_path.unlink()
                        return False
                    f.write(chunk)

            self.state.downloaded_files.add(str(output_path))
            self.log(f"Downloaded: {output_path} ({size} bytes)")
            return True

        except requests.RequestException as e:
            self.log(f"Failed to download {url}: {e}", "ERROR")
            return False

    def html_to_markdown(self, html: str, url: str) -> str:
        """Convert HTML to markdown, preserving tables"""
        soup = BeautifulSoup(html, "html.parser")

        # Remove script, style, and navigation elements
        for element in soup(["script", "style", "noscript"]):
            element.decompose()

        # Remove tree menu (navigation sidebar)
        for element in soup.find_all("div", class_="dtree"):
            element.decompose()
        for element in soup.find_all("div", id="treemenu"):
            element.decompose()

        lines = []
        lines.append(f"<!-- Source: {url} -->")
        lines.append(f"<!-- Crawled: {time.strftime('%Y-%m-%d %H:%M:%S')} -->")
        lines.append("")

        # Get title
        title = soup.find("title")
        if title:
            title_text = title.get_text().strip()
            # Clean up title
            title_text = re.sub(r'\s+', ' ', title_text)
            lines.append(f"# {title_text}")
            lines.append("")

        # For ligand results pages, extract the specific data table
        if "ligand_results.cgi" in url:
            lines.extend(self._extract_ligand_results(soup, url))
        else:
            # Find main content area - look for specific CASP content divs
            content = (
                soup.find("div", id="content") or
                soup.find("div", class_="content") or
                soup.find("td", class_="content") or
                soup.find("main") or
                soup.body
            )
            if not content:
                content = soup

            # Process data tables (not navigation tables)
            for table in content.find_all("table"):
                # Skip small navigation tables
                rows = table.find_all("tr")
                if len(rows) > 2:  # Only process tables with actual data
                    table_md = self._table_to_markdown(table)
                    if table_md.strip():
                        lines.append(table_md)
                        lines.append("")

            # Get remaining text content
            text = content.get_text(separator="\n", strip=True)
            # Clean up excessive whitespace
            text = re.sub(r"\n{3,}", "\n\n", text)
            # Remove duplicate lines
            seen = set()
            clean_lines = []
            for line in text.split("\n"):
                line = line.strip()
                if line and line not in seen and len(line) > 3:
                    seen.add(line)
                    clean_lines.append(line)
            lines.append("\n".join(clean_lines))

        return "\n".join(lines)

    def _extract_ligand_results(self, soup, url: str) -> list:
        """Extract ligand results data specifically"""
        lines = []

        # Find target name
        target_match = re.search(r'target=(\w+)', url)
        if target_match:
            lines.append(f"## Target: {target_match.group(1)}")
            lines.append("")

        # Find the results table - it has LDDT_pli, RMSD columns
        for table in soup.find_all("table"):
            headers = [th.get_text().strip() for th in table.find_all("th")]
            # Check if this is the results table
            if any("LDDT" in h or "RMSD" in h for h in headers):
                lines.append("### Prediction Results")
                lines.append("")
                lines.append(self._table_to_markdown(table))
                lines.append("")
                break

        # Also look for tables with class or id suggesting results
        for table in soup.find_all("table", class_=re.compile(r"result|data|score", re.I)):
            table_md = self._table_to_markdown(table)
            if table_md.strip() and table_md not in "\n".join(lines):
                lines.append(table_md)
                lines.append("")

        return lines

    def _table_to_markdown(self, table) -> str:
        """Convert HTML table to markdown table"""
        rows = []
        headers = []

        # Get headers from th elements
        header_row = table.find("tr")
        if header_row:
            for th in header_row.find_all("th"):
                text = th.get_text().strip()
                text = re.sub(r'\s+', ' ', text)  # Normalize whitespace
                if text:
                    headers.append(text)

        # If no th elements, try first row td elements
        if not headers and header_row:
            for td in header_row.find_all("td"):
                text = td.get_text().strip()
                text = re.sub(r'\s+', ' ', text)
                if text:
                    headers.append(text)

        if headers:
            # Clean headers - remove very long repeated text
            clean_headers = []
            for h in headers:
                if len(h) > 50:
                    h = h[:47] + "..."
                clean_headers.append(h)
            headers = clean_headers

            rows.append("| " + " | ".join(headers) + " |")
            rows.append("| " + " | ".join(["---"] * len(headers)) + " |")

        # Get data rows
        data_rows = table.find_all("tr")[1:] if headers else table.find_all("tr")
        for tr in data_rows:
            cells = []
            for td in tr.find_all("td"):
                text = td.get_text().strip()
                text = re.sub(r'\s+', ' ', text)  # Normalize whitespace
                text = text.replace("|", "\\|")
                # Truncate very long cells
                if len(text) > 50:
                    text = text[:47] + "..."
                cells.append(text)

            # Skip empty rows or rows with all empty cells
            if cells and any(c.strip() for c in cells):
                # Pad cells to match header count
                if headers:
                    while len(cells) < len(headers):
                        cells.append("")
                    cells = cells[:len(headers)]
                rows.append("| " + " | ".join(cells) + " |")

        return "\n".join(rows)

    def save_page(self, url: str, content: str, output_path: Path):
        """Save page content as markdown"""
        if self.dry_run:
            self.log(f"Would save: {output_path}")
            return

        output_path.parent.mkdir(parents=True, exist_ok=True)
        md = self.html_to_markdown(content, url)
        with open(output_path, "w") as f:
            f.write(md)
        self.log(f"Saved: {output_path}")

    def crawl_casp16_ligand(self):
        """Crawl CASP16 ligand results"""
        self.log("=== Crawling CASP16 Ligand Results ===")

        # Main ligand results page
        url = f"{BASE_URL}/casp16/results.cgi?tr_type=ligand"
        resp = self.fetch(url)
        if resp:
            output = OUTPUT_DIR / "casp16" / "results" / "ligand" / "index.md"
            self.save_page(url, resp.text, output)

        # Individual target results
        for target in CASP16_LIGAND_TARGETS:
            url = f"{BASE_URL}/casp16/ligand_results.cgi?target={target}"
            resp = self.fetch(url)
            if resp:
                output = OUTPUT_DIR / "casp16" / "results" / "ligand" / f"{target}.md"
                self.save_page(url, resp.text, output)

        # Z-scores page
        url = f"{BASE_URL}/casp16/zscores_ligand.cgi"
        resp = self.fetch(url)
        if resp:
            output = OUTPUT_DIR / "casp16" / "rankings" / "ligand.md"
            self.save_page(url, resp.text, output)

    def crawl_casp15_ligand(self):
        """Crawl CASP15 ligand results"""
        self.log("=== Crawling CASP15 Ligand Results ===")

        # Main ligand results page
        url = f"{BASE_URL}/casp15/results.cgi?tr_type=ligand"
        resp = self.fetch(url)
        if resp:
            output = OUTPUT_DIR / "casp15" / "results" / "ligand" / "index.md"
            self.save_page(url, resp.text, output)

        # Individual target results
        for target in CASP15_LIGAND_TARGETS:
            url = f"{BASE_URL}/casp15/ligand_results.cgi?target={target}"
            resp = self.fetch(url)
            if resp:
                output = OUTPUT_DIR / "casp15" / "results" / "ligand" / f"{target}.md"
                self.save_page(url, resp.text, output)

        # Z-scores page
        url = f"{BASE_URL}/casp15/zscores_ligand.cgi"
        resp = self.fetch(url)
        if resp:
            output = OUTPUT_DIR / "casp15" / "rankings" / "ligand.md"
            self.save_page(url, resp.text, output)

    def crawl_casp16_overview(self):
        """Crawl CASP16 overview pages"""
        self.log("=== Crawling CASP16 Overview ===")

        pages = [
            ("/casp16/", "casp16/index.md"),
            ("/casp16/results.cgi", "casp16/results/overview.md"),
        ]

        for path, output in pages:
            url = f"{BASE_URL}{path}"
            resp = self.fetch(url)
            if resp:
                self.save_page(url, resp.text, OUTPUT_DIR / output)

    def crawl_casp15_overview(self):
        """Crawl CASP15 overview pages"""
        self.log("=== Crawling CASP15 Overview ===")

        pages = [
            ("/casp15/", "casp15/index.md"),
            ("/casp15/results.cgi", "casp15/results/overview.md"),
        ]

        for path, output in pages:
            url = f"{BASE_URL}{path}"
            resp = self.fetch(url)
            if resp:
                self.save_page(url, resp.text, OUTPUT_DIR / output)

    def download_sequences(self):
        """Download sequence files"""
        self.log("=== Downloading Sequence Files ===")

        for filename in CASP16_SEQUENCE_FILES:
            url = f"{BASE_URL}/download_area/CASP16/sequences/{filename}"
            output = OUTPUT_DIR / "casp16" / "sequences" / filename
            self.download_file(url, output)

    def download_pharma_ligands(self):
        """Download pharma ligand metadata files"""
        self.log("=== Downloading Pharma Ligand Files ===")

        for filename in PHARMA_LIGAND_FILES:
            url = f"{BASE_URL}/download_area/CASP16/targets/pharma_ligands/{filename}"
            output = OUTPUT_DIR / "casp16" / "pharma_ligands" / filename
            # Use smaller limit for tar.gz files that might be larger
            max_size = 100 * 1024 if filename.endswith(".tar.gz") else MAX_FILE_SIZE
            self.download_file(url, output, max_size=max_size)

    def download_ligand_results(self):
        """Download ligand scoring/results CSV files"""
        self.log("=== Downloading Ligand Scoring Results ===")

        for remote_path, local_path, max_size in LIGAND_RESULTS_FILES:
            url = f"{BASE_URL}/download_area/{remote_path}"
            output = OUTPUT_DIR / local_path
            self.download_file(url, output, max_size=max_size)

    def download_ligand_predictions(self, include_casp15: bool = True):
        """Download ligand prediction archives"""
        self.log("=== Downloading CASP16 Ligand Predictions ===")

        for filename, max_size in CASP16_LIGAND_PREDICTIONS:
            url = f"{BASE_URL}/download_area/CASP16/predictions/ligands/{filename}"
            output = OUTPUT_DIR / "casp16" / "predictions" / "ligands" / filename
            self.download_file(url, output, max_size=max_size)

        # Pharma affinity files (small)
        self.log("=== Downloading CASP16 Pharma Affinity Files ===")
        for filename, max_size in CASP16_PHARMA_AFFINITY:
            url = f"{BASE_URL}/download_area/CASP16/predictions/ligands/pharma/{filename}"
            output = OUTPUT_DIR / "casp16" / "predictions" / "pharma" / filename
            self.download_file(url, output, max_size=max_size)

        if include_casp15:
            self.log("=== Downloading CASP15 Ligand Predictions ===")
            for filename, max_size in CASP15_LIGAND_PREDICTIONS:
                url = f"{BASE_URL}/download_area/CASP15/predictions/ligands/{filename}"
                output = OUTPUT_DIR / "casp15" / "predictions" / "ligands" / filename
                self.download_file(url, output, max_size=max_size)

    def run(self, section: Optional[str] = None, include_casp15_predictions: bool = True):
        """Run the spider"""
        self.log(f"Starting CASP Spider (dry_run={self.dry_run})")

        try:
            if section is None or section == "ligand":
                self.crawl_casp16_ligand()
                self.crawl_casp15_ligand()

            if section is None or section == "overview":
                self.crawl_casp16_overview()
                self.crawl_casp15_overview()

            if section is None or section == "downloads":
                self.download_sequences()
                self.download_pharma_ligands()

            if section is None or section == "results":
                self.download_ligand_results()

            if section is None or section == "predictions":
                self.download_ligand_predictions(include_casp15=include_casp15_predictions)

            # Save state
            if not self.dry_run:
                self.state.save(STATE_FILE)
                self.log(f"State saved to {STATE_FILE}")

            self.log("=== Spider Complete ===")
            self.log(f"Crawled {len(self.state.crawled_urls)} URLs")
            self.log(f"Downloaded {len(self.state.downloaded_files)} files")
            if self.state.failed_urls:
                self.log(f"Failed: {len(self.state.failed_urls)} URLs", "WARN")

        except KeyboardInterrupt:
            self.log("Interrupted by user", "WARN")
            if not self.dry_run:
                self.state.save(STATE_FILE)
                self.log(f"State saved - use --resume to continue")


def main():
    parser = argparse.ArgumentParser(description="Spider CASP15/16 ligand prediction data")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be crawled without fetching")
    parser.add_argument("--delay", type=float, default=DEFAULT_DELAY, help="Delay between requests (seconds)")
    parser.add_argument("--section", choices=["ligand", "overview", "downloads", "results", "predictions"],
                        help="Only crawl specific section")
    parser.add_argument("--resume", action="store_true", help="Resume from saved state")
    parser.add_argument("--skip-casp15-predictions", action="store_true",
                        help="Skip downloading CASP15 prediction archives (saves ~450MB)")

    args = parser.parse_args()

    spider = CASPSpider(
        dry_run=args.dry_run,
        delay=args.delay,
        resume=args.resume,
    )
    spider.run(
        section=args.section,
        include_casp15_predictions=not args.skip_casp15_predictions
    )


if __name__ == "__main__":
    main()

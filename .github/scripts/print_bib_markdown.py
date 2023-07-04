import os
import requests
import logging
from typing import Any, Optional


BASE_URL = "https://api.adsabs.harvard.edu/v1"
ALPS_ADS_BIB_CODE = "2018JPlPh..84d9003V"
ALPS_LIBRARY_ID = "rrfGfF0OT4CH7y4AOOil4w"   # TODO: change

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d %(levelname)s - %(funcName)s: %(message)s',
    datefmt='%H:%M:%S',
)
logger = logging.getLogger(__name__)


def data_from_library_api(**params: Any) -> dict:
    logger.info(f"Getting library data for libray (id: {ALPS_LIBRARY_ID}) with {params}")

    response = requests.get(
        url=f"{BASE_URL}/biblib/libraries/{ALPS_LIBRARY_ID}",
        headers={'Authorization': f"Bearer {os.environ['ADS_API_TOKEN']}"},
        params=params
    )

    if response.status_code != 200:
        raise RuntimeError(f"Failed to get {response.url} due to:\n {response.text}")

    return response.json()


def get_bibcodes_in_library() -> list:
    logger.info("Getting a list of bibcodes in an ADS library")

    data = data_from_library_api(rows=0)
    num_papers_in_library = data['metadata']['num_documents']
    logger.info(f"Found {num_papers_in_library} papers in library")

    data = data_from_library_api(rows=10000, fl="bibcode")
    return data["documents"]


def data_from_export_api(**params: Any) -> dict:
    response = requests.post(
        url=f"{BASE_URL}/export/custom",
        headers={'Authorization': f"Bearer {os.environ['ADS_API_TOKEN']}"},
        json=params
    )

    if response.status_code != 200:
        raise RuntimeError(f"Failed to post {response.url} due to:\n {response.text}")

    return response.json()


def entry_from_bibcode(bibcode: str) -> Optional[dict]:

    # See https://adsabs.github.io/help/actions/export for formatting options
    seperator = ":|:"
    string_format = seperator.join(["%3.2l", "%Y", "%J", "%d", "%T"])
    data = data_from_export_api(bibcode=[bibcode], format=string_format)

    try:
        authors, year, journal, doi, title = data["export"].split(seperator)
    except IndexError as e:
        logger.error(f"Failed to parse: {data}.\nException: {e}")
        return None

    entry = {
        'bibcode': bibcode,
        'authors': authors,
        'title': title.strip("\n"),
        'year': int(year),
        'doi': doi,
        'doiurl': f"https://dx.doi.org/{doi}",
        'adsurl': f"https://ui.adsabs.harvard.edu/abs/{bibcode}",
        'journal': journal
    }
    return entry


def get_sorted_bib_entries() -> list:
    logger.info("Getting sorted bib entries by year (high->low) and author list (a->z)")

    entries = [entry_from_bibcode(c) for c in get_bibcodes_in_library()]
    return sorted(
        [entry for entry in entries if entry is not None],
        key=lambda x: (-x["year"], x["authors"])
    )


def print_markdown() -> None:
    year: Optional[int] = None

    print(
        "title: Papers citing ALPS\n"
        "# Papers citing ALPS"
    )

    for entry in get_sorted_bib_entries():
        if entry["year"] != year:
            year = entry["year"]
            print(f"## {year}")

        print(f"{entry['authors']}, _{entry['title']}_, {entry['journal']}, [{entry['doi']}]({entry['doiurl']})\n")


def main() -> None:

    if "ADS_API_TOKEN" not in os.environ:
        raise RuntimeError("Please set ADS_API_TOKEN as an environment variable")

    print_markdown()


if __name__ == '__main__':
    main()

import os
import json
import requests

BASE_URL = "https://api.adsabs.harvard.edu/v1"
ADS_LIBRARY_ID = "vPWlMyiDRq-YfZ2bGDWS7Q"

BIBCODES = [
    "2018JPlPh..84d9003V", #2018 JPP
    "2025PhPl...32i2104K"  #2025 PoP
]

def query_response() -> requests.Response:
    """Get a response by querying the ADS database for citations of the ALPS code"""

    # build the OR-joined clause:
    or_clause = " OR ".join(BIBCODES)
    query = f"citations(bibcode:({or_clause}))"
    
    response = requests.get(
        url=f"{BASE_URL}/search/query",
        headers={"Authorization": f"Bearer {os.environ['ADS_API_TOKEN']}"},
        params={
            "q": query,  # Query
            "fl": "bibcode",  # Fields to return
            "rows": 1000  # Max num results to return
        }
    )
    return response


def bibcodes_from_response(response: requests.Response) -> list:
    """Given a requests response extract the bibcode unique identifiers"""

    if response.status_code != 200:
        raise RuntimeError(f"Failed to get: {response.url} content:\n{response.content}")

    data = json.loads(response.content.decode())

    try:
        bibcodes = [item["bibcode"] for item in data["response"]["docs"]]
        assert len(bibcodes) == data["response"]["numFound"]
        print(f"Found {len(bibcodes)} papers citing {ALPS_ADS_BIB_CODE}")

    except (KeyError, AssertionError) as e:
        raise RuntimeError(f"Response from {response.url} was not valid") from e

    return bibcodes


def add_bibcodes_to_library(bibcodes: list, library_id: str) -> None:
    """Add a list of papers to an ADS library"""

    response = requests.post(
        url=f"{BASE_URL}/biblib/documents/{library_id}",
        headers={"Authorization": f"Bearer {os.environ['ADS_API_TOKEN']}"},
        json={"bibcode": bibcodes, "action": "add"}
    )

    if response.status_code != 200:
        raise RuntimeError(f"Failed to post {response.url} content:\n {response.content}")

    data = json.loads(response.content.decode())
    print(f"Added {data['number_added']} papers to library (id: {library_id})")
    return None


def main() -> None:

    if "ADS_API_TOKEN" not in os.environ:
        raise RuntimeError("Please set ADS_API_TOKEN as an environment variable")

    add_bibcodes_to_library(
        bibcodes=bibcodes_from_response(query_response()),
        library_id=ADS_LIBRARY_ID
    )


if __name__ == '__main__':
    main()

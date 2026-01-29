import requests
import urllib.parse

def search_ols(term, ontology=None, rows=1):
    base_url = "https://www.ebi.ac.uk/ols/api/search"
    params = {
        "q": term,
        "rows": rows,
        "exact": "false",
        "type": "class"
    }
    if ontology:
        params["ontology"] = ontology

    print(f"Querying: {term} (ontology={ontology})...")
    response = requests.get(base_url, params=params)
    
    if response.status_code != 200:
        print(f"Error: {response.status_code}")
        return

    data = response.json()
    docs = data.get("response", {}).get("docs", [])

    if not docs:
        print("No results found.")
        return

    for doc in docs:
        print(f"  Match: {doc.get('label')} ({doc.get('obo_id')})")
        print(f"  Ontology: {doc.get('ontology_name')}")
        print(f"  Score: {doc.get('score')}")
        print("-" * 20)

if __name__ == "__main__":
    print("Testing UBERON (Tissue)...")
    search_ols("liver", "uberon")
    
    print("\nTesting CL (Cell Type)...")
    search_ols("neuron", "cl")
    
    print("\nTesting NCBITaxon (Organism)...")
    search_ols("homo sapiens", "ncbitaxon")
    
    print("\nTesting Fuzzy...")
    search_ols("lung tissue", "uberon")

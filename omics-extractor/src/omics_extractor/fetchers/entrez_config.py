"""Centralized Entrez configuration and polite rate limiting.

Reads NCBI settings from environment variables:
- NCBI_EMAIL: required by NCBI; empty string if not set
- NCBI_API_KEY: optional E-utilities key for higher limits
- NCBI_TOOL: tool name reported to NCBI (default: omics-extractor)
- NCBI_SLEEP: seconds to sleep between requests (default: 0.34)
"""

from Bio import Entrez
import os
import time
import random
from urllib.error import HTTPError

EMAIL = os.getenv('NCBI_EMAIL', '')
API_KEY = os.getenv('NCBI_API_KEY')
TOOL = os.getenv('NCBI_TOOL', 'omics-extractor')
_SLEEP = float(os.getenv('NCBI_SLEEP', '0.34'))

_configured = False

def configure() -> None:
    """Apply Entrez global configuration once per process."""
    global _configured
    if _configured:
        return
    Entrez.email = EMAIL
    Entrez.tool = TOOL
    if API_KEY:
        Entrez.api_key = API_KEY
    _configured = True


def rate_limit() -> None:
    """Polite delay to avoid hammering NCBI endpoints."""
    if _SLEEP > 0:
        import time as _t
        _t.sleep(_SLEEP)


def _retryable(code: int) -> bool:
    return code in (429, 500, 502, 503, 504)


def backoff_sleep(attempt: int) -> None:
    """Exponential backoff with jitter, capped."""
    base = max(_SLEEP, 0.2)
    delay = min(10.0, base * (2 ** attempt) + random.uniform(0, 0.3))
    time.sleep(delay)


def with_retries(call, max_retries: int = 5, read: bool = False):
    """Call an Entrez function with retry/backoff on 429/5xx.

    If read=True, parse the handle via Entrez.read() and return the parsed object.
    """
    last_exc = None
    for attempt in range(max_retries):
        try:
            handle = call()
            if read:
                out = Entrez.read(handle)
                try:
                    handle.close()
                except Exception:
                    pass
                rate_limit()
                return out
            else:
                data = handle.read()
                try:
                    handle.close()
                except Exception:
                    pass
                rate_limit()
                return data
        except HTTPError as e:
            last_exc = e
            if _retryable(getattr(e, 'code', 0)):
                backoff_sleep(attempt)
                continue
            raise
        except Exception as e:
            last_exc = e
            backoff_sleep(attempt)
            continue
    if last_exc:
        raise last_exc
    raise RuntimeError("Entrez request failed without exception")

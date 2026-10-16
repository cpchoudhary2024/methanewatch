# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Authenticated granule download from CDSE, cached under $METHANE_CACHE_DIR."""

from __future__ import annotations

import os
import time
from pathlib import Path

import requests

from .catalog import Granule

_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
)
_DOWNLOAD_URL = "https://download.dataspace.copernicus.eu/odata/v1/Products({pid})/$value"
_TIMEOUT = 60


class CredentialsError(RuntimeError):
    pass


def _load_dotenv(path: str = ".env") -> None:
    p = Path(path)
    if not p.exists():
        return
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        os.environ.setdefault(k.strip(), v.strip())


_token: dict = {"value": None, "expires": 0.0}


def get_access_token(username: str | None = None, password: str | None = None) -> str:
    """A CDSE access token, reused until a minute before it expires."""
    if _token["value"] and time.time() < _token["expires"] - 60:
        return _token["value"]
    _load_dotenv()
    username = username or os.environ.get("CDSE_USERNAME", "")
    password = password or os.environ.get("CDSE_PASSWORD", "")
    if not username or not password:
        raise CredentialsError(
            "CDSE_USERNAME / CDSE_PASSWORD not set. Register at "
            "https://dataspace.copernicus.eu/ and put them in .env (see .env.example)."
        )
    resp = requests.post(
        _TOKEN_URL,
        data={
            "grant_type": "password",
            "username": username,
            "password": password,
            "client_id": "cdse-public",
        },
        timeout=_TIMEOUT,
    )
    if resp.status_code == 401:
        raise CredentialsError(
            "CDSE rejected the credentials (401). Double-check CDSE_USERNAME / "
            "CDSE_PASSWORD and that you confirmed the account email."
        )
    resp.raise_for_status()
    body = resp.json()
    _token.update(value=body["access_token"], expires=time.time() + body.get("expires_in", 600))
    return _token["value"]


def cache_dir() -> Path:
    d = Path(os.environ.get("METHANE_CACHE_DIR", "data/cache"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def download_granule(
    granule: Granule,
    token: str | None = None,
    *,
    force: bool = False,
) -> Path:
    dest = cache_dir() / granule.name
    if dest.exists() and not force and dest.stat().st_size > 0:
        return dest

    if not granule.online:
        raise RuntimeError(
            f"{granule.name} is in CDSE long-term archive (offline). "
            "Trigger an LTA retrieval and retry later, or pick another date."
        )

    token = token or get_access_token()
    url = _DOWNLOAD_URL.format(pid=granule.product_id)
    with requests.get(
        url,
        headers={"Authorization": f"Bearer {token}"},
        stream=True,
        timeout=_TIMEOUT,
        allow_redirects=True,
    ) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
        tmp.rename(dest)
    return dest

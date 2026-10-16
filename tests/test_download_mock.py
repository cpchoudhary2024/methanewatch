from datetime import datetime, timezone
from unittest.mock import MagicMock

import pytest
import requests

from methane.catalog import Granule
from methane.download import (
    CredentialsError,
    _load_dotenv,
    cache_dir,
    download_granule,
    get_access_token,
)


@pytest.fixture(autouse=True)
def _fresh_token():
    from methane import download

    download._token.update(value=None, expires=0.0)


def _granule(name="S5P_OFFL_L2__CH4____TEST.nc", online=True, pid="pid-1"):
    return Granule(
        product_id=pid,
        name=name,
        start=datetime(2024, 6, 1, tzinfo=timezone.utc),
        size_bytes=1024,
        online=online,
    )


def test_get_access_token_raises_without_credentials(monkeypatch, tmp_path, mocker):
    monkeypatch.delenv("CDSE_USERNAME", raising=False)
    monkeypatch.delenv("CDSE_PASSWORD", raising=False)
    monkeypatch.chdir(tmp_path)  # no .env here either
    mock_post = mocker.patch("requests.post")

    with pytest.raises(CredentialsError, match="CDSE_USERNAME"):
        get_access_token()
    mock_post.assert_not_called()


def test_get_access_token_success(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"access_token": "fake-token-abc"}
    mock_resp.raise_for_status.return_value = None
    mock_post = mocker.patch("requests.post", return_value=mock_resp)

    token = get_access_token(username="u@example.com", password="pw")
    assert token == "fake-token-abc"
    mock_post.assert_called_once()
    assert mock_post.call_args.kwargs["data"]["grant_type"] == "password"


def test_get_access_token_401_raises_credentials_error(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 401
    mocker.patch("requests.post", return_value=mock_resp)

    with pytest.raises(CredentialsError, match="rejected the credentials"):
        get_access_token(username="u@example.com", password="wrong")


def test_get_access_token_sudden_network_dropout(mocker):
    """The CDSE identity server connection drops mid-request."""
    mocker.patch("requests.post", side_effect=requests.exceptions.ConnectionError("dropped"))
    with pytest.raises(requests.exceptions.ConnectionError):
        get_access_token(username="u@example.com", password="pw")


def test_load_dotenv_does_not_override_real_env(tmp_path, monkeypatch):
    envfile = tmp_path / ".env"
    envfile.write_text("CDSE_USERNAME=fromfile\nCDSE_PASSWORD=alsofromfile\n# comment\n\n")
    monkeypatch.setenv("CDSE_USERNAME", "already-set")
    monkeypatch.delenv("CDSE_PASSWORD", raising=False)

    _load_dotenv(str(envfile))

    import os

    assert os.environ["CDSE_USERNAME"] == "already-set"  # real env wins
    assert os.environ["CDSE_PASSWORD"] == "alsofromfile"  # .env fills the gap


def test_load_dotenv_missing_file_is_a_noop(tmp_path):
    # Must not raise just because .env doesn't exist (e.g. CI/Docker with
    # credentials injected purely via env vars).
    _load_dotenv(str(tmp_path / "does_not_exist.env"))


def test_cache_dir_creates_directory(tmp_path, monkeypatch):
    target = tmp_path / "nested" / "cache"
    monkeypatch.setenv("METHANE_CACHE_DIR", str(target))
    d = cache_dir()
    assert d == target
    assert d.is_dir()


def test_download_granule_is_idempotent_on_existing_cache(tmp_path, monkeypatch, mocker):
    monkeypatch.setenv("METHANE_CACHE_DIR", str(tmp_path))
    g = _granule()
    existing = tmp_path / g.name
    existing.write_bytes(b"already here")

    mock_get = mocker.patch("requests.get")
    path = download_granule(g, token="tok")

    assert path == existing
    mock_get.assert_not_called()  # cache hit -> no network call at all


def test_download_granule_rejects_offline_granule(tmp_path, monkeypatch):
    monkeypatch.setenv("METHANE_CACHE_DIR", str(tmp_path))
    g = _granule(online=False)
    with pytest.raises(RuntimeError, match="long-term archive"):
        download_granule(g, token="tok")


def test_download_granule_streams_to_cache(tmp_path, monkeypatch, mocker):
    monkeypatch.setenv("METHANE_CACHE_DIR", str(tmp_path))
    g = _granule(name="fresh_granule.nc")

    chunks = [b"chunk-one-", b"chunk-two-", b"chunk-three"]
    mock_stream = MagicMock()
    mock_stream.raise_for_status.return_value = None
    mock_stream.iter_content.return_value = iter(chunks)
    mock_ctx = MagicMock()
    mock_ctx.__enter__.return_value = mock_stream
    mock_ctx.__exit__.return_value = False
    mocker.patch("requests.get", return_value=mock_ctx)

    path = download_granule(g, token="tok")
    assert path.exists()
    assert path.read_bytes() == b"".join(chunks)
    assert not path.with_suffix(path.suffix + ".part").exists()  # temp file renamed away


def test_download_granule_sudden_network_dropout_leaves_no_partial_final_file(
    tmp_path, monkeypatch, mocker
):
    """A connection reset mid-stream must not leave a corrupt file at the
    final destination path, only the .part temp file, if anything."""
    monkeypatch.setenv("METHANE_CACHE_DIR", str(tmp_path))
    g = _granule(name="interrupted.nc")

    mock_stream = MagicMock()
    mock_stream.raise_for_status.return_value = None

    def _boom(*a, **k):
        yield b"partial-bytes"
        raise requests.exceptions.ConnectionError("connection reset")

    mock_stream.iter_content.side_effect = _boom
    mock_ctx = MagicMock()
    mock_ctx.__enter__.return_value = mock_stream
    mock_ctx.__exit__.return_value = False
    mocker.patch("requests.get", return_value=mock_ctx)

    with pytest.raises(requests.exceptions.ConnectionError):
        download_granule(g, token="tok")

    dest = tmp_path / g.name
    assert not dest.exists()  # never renamed into place


def test_download_granule_fetches_token_when_not_supplied(tmp_path, monkeypatch, mocker):
    monkeypatch.setenv("METHANE_CACHE_DIR", str(tmp_path))
    g = _granule(name="needs_token.nc")
    mock_get_token = mocker.patch("methane.download.get_access_token", return_value="auto-token")

    mock_stream = MagicMock()
    mock_stream.raise_for_status.return_value = None
    mock_stream.iter_content.return_value = iter([b"data"])
    mock_ctx = MagicMock()
    mock_ctx.__enter__.return_value = mock_stream
    mock_ctx.__exit__.return_value = False
    mock_get = mocker.patch("requests.get", return_value=mock_ctx)

    download_granule(g, token=None)
    mock_get_token.assert_called_once()
    assert mock_get.call_args.kwargs["headers"]["Authorization"] == "Bearer auto-token"

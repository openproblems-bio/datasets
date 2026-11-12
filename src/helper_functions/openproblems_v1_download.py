import shutil
import urllib.request


def download_url(url, destination):
  """Download a file from a URL, without a spoofed User-Agent.

  Drop-in replacement for `scprep.io.download.download_url`, which sends
  `User-Agent: Mozilla/5.0` and installs it as the global urllib opener.
  figshare and GEO refuse such requests nowadays.
  """
  if isinstance(destination, str):
    with open(destination, "wb") as handle:
      return download_url(url, handle)
  with urllib.request.urlopen(url) as response:
    if response.status != 200:
      raise RuntimeError(f"Downloading {url} returned HTTP {response.status} instead of the file")
    shutil.copyfileobj(response, destination)


def use_plain_downloads():
  """Make openproblems v0.8.0 download its datasets with `download_url()`."""
  import scprep.io.download

  scprep.io.download.download_url = download_url

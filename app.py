"""mmcifXplorer — run with:  python app.py  [--port 8050] [--debug]"""
from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore", message="You may be importing Biopython")

import dash  # noqa: E402

from mmcifxplorer.callbacks import register  # noqa: E402
from mmcifxplorer.layout import serve_layout  # noqa: E402

app = dash.Dash(__name__, title="mmcifXplorer", suppress_callback_exceptions=True,
                update_title=None)
app.layout = serve_layout
register(app)
server = app.server

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="mmcifXplorer dashboard")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8050)
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args()
    app.run(host=a.host, port=a.port, debug=a.debug)

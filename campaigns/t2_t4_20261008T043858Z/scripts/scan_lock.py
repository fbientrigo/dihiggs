"""Nonblocking interprocess lock for mutating campaign scan output."""
from __future__ import annotations

import fcntl
import json
import os
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def campaign_scan_lock(path: Path, identity: str):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    handle=path.open("a+")
    locked=False
    try:
        try:
            fcntl.flock(handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner=handle.read().strip() or "owner metadata unavailable"
            raise RuntimeError(f"T4 scan lock is held: {owner}") from exc
        locked=True
        handle.seek(0)
        handle.truncate()
        handle.write(json.dumps({"pid":os.getpid(),"scan_identity":identity},sort_keys=True)+"\n")
        handle.flush()
        os.fsync(handle.fileno())
        yield handle
    finally:
        if not handle.closed:
            if locked:
                try:
                    handle.seek(0)
                    handle.truncate()
                    handle.flush()
                    fcntl.flock(handle.fileno(),fcntl.LOCK_UN)
                finally:
                    handle.close()
            else:
                handle.close()

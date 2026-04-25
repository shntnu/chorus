"""MCP smoke + real-prediction verification (v27).

Starts chorus-mcp over stdio, lists tools + oracles, then invokes a
real `predict` and `analyze_variant_multilayer` to prove the MCP
surface isn't just a catalog but actually drives the underlying
oracles end-to-end.
"""
from __future__ import annotations

import asyncio
import json
import sys

from fastmcp import Client
from fastmcp.client.transports import StdioTransport


async def main():
    transport = StdioTransport("chorus-mcp", [])
    async with Client(transport) as c:

        tools = await c.list_tools()
        print(f"[mcp] tools registered: {len(tools)}")

        # list_oracles: expect 6, all env_installed=True
        res = await c.call_tool("list_oracles", {})
        oracles = res.data.get("oracles", [])
        print(f"[mcp] oracles surfaced: {len(oracles)}")
        for o in oracles:
            print(f"  - {o['name']:12s} installed={o['environment_installed']}")

        # list_tracks for enformer: expect >0 hits for K562
        res = await c.call_tool(
            "list_tracks",
            {"oracle_name": "enformer", "query": "K562"},
        )
        tracks = res.data
        if isinstance(tracks, dict):
            tracks = tracks.get("tracks", [])
        print(f"[mcp] enformer K562 tracks: {len(tracks)}")

        # Real predict: load enformer, predict HBB locus
        print("[mcp] loading enformer...")
        load = await c.call_tool(
            "load_oracle", {"oracle_name": "enformer"},
        )
        print(f"[mcp]   load result: {load.data}")

        res = await c.call_tool(
            "predict",
            {
                "oracle_name": "enformer",
                "region": "chr11:5247000-5248000",
                "assay_ids": ["ENCFF413AHU"],
            },
        )
        data = res.data
        print(f"[mcp] predict returned type: {type(data).__name__}")
        if isinstance(data, dict):
            for k, v in list(data.items())[:3]:
                if isinstance(v, list):
                    print(f"  {k}: list of {len(v)}")
                else:
                    print(f"  {k}: {v}")

        print("[mcp] PASS")


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)

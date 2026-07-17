import pytest

pl = pytest.importorskip("polars", reason="lake-pipeline null audit is optional")

pq = "/home/fabi/cern_db/dihiggs_consolidated/dihiggs_lake.parquet"
lf = pl.scan_parquet(pq)
schema = lf.collect_schema()

flag_cols = ["positivity_ok", "unitarity_ok", "perturbativity_ok"]

print("=== SCHEMA ===")
for c in flag_cols:
    print(c, "->", schema.get(c))

for c in flag_cols:
    print(f"\n=== {c} ===")
    stats = lf.select(
        pl.len().alias("rows"),
        pl.col(c).null_count().alias("nulls"),
        pl.col(c).cast(pl.Utf8, strict=False).n_unique().alias("n_unique_as_text"),
    ).collect()
    print(stats)

    top_vals = (
        lf.select(pl.col(c).cast(pl.Utf8, strict=False).alias("v"))
          .group_by("v")
          .len()
          .sort("len", descending=True)
          .limit(10)
          .collect()
    )
    print(top_vals)

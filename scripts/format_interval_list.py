#!/usr/bin/python


with open("master_interval_list.list", "w") as out:
    with open("picard_interval_list.list", "r") as inp:
        for line in inp:
            if not line.startswith("@"):
                line = line.strip().split("\t")
                chrom, start, end = line[0], line[1], line[2]
                print(f"{chrom}:{start}-{end}", file=out)
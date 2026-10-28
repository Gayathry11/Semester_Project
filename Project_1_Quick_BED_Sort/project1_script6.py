#!/usr/bin/env python3

import sys


def read_chromosome_order(selection_file):
    """
    Read chromosome names from the selection file and
    preserve their order.
    """

    chromosome_order = []

    with open(selection_file, "r") as infile:
        for line in infile:

            chromosome = line.strip()

            if chromosome:
                chromosome_order.append(chromosome)

    return chromosome_order


def initialize_storage(chromosome_order):
    """
    Create an empty list for each chromosome.
    """

    return {chromosome: [] for chromosome in chromosome_order}


def store_bed_records(chromosome_storage):
    """
    Read BED records from standard input and
    store them chromosome-wise.
    """

    for line in sys.stdin:

        cleaned_line = line.strip()

        if not cleaned_line:
            continue

        columns = cleaned_line.split("\t")

        if len(columns) < 3:
            continue

        chromosome_name = columns[0]

        if chromosome_name not in chromosome_storage:
            continue

        try:
            start_coordinate = int(columns[1])
        except ValueError:
            continue

        chromosome_storage[chromosome_name].append(
            (start_coordinate, line)
        )


def write_sorted_records(chromosome_order, chromosome_storage):
    """
    Output chromosomes in selection-file order and
    sort entries by genomic start coordinate.
    """

    for chromosome in chromosome_order:

        records = chromosome_storage[chromosome]

        records.sort(key=lambda record: record[0])

        for _, bed_record in records:
            sys.stdout.write(bed_record)


def main():

    selection_file = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "standard_selection.tsv"
    )

    chromosome_order = read_chromosome_order(
        selection_file
    )

    chromosome_storage = initialize_storage(
        chromosome_order
    )

    store_bed_records(chromosome_storage)

    write_sorted_records(
        chromosome_order,
        chromosome_storage
    )


if __name__ == "__main__":
    main()

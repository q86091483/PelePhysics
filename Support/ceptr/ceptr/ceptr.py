"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib

import cantera as ct

import ceptr.converter as converter


def convert(
    fname,
    hformat,
    remove_1,
    remove_pow,
    min_op_count,
    recursive_op_count,
    store_in_jacobian,
    round_decimals,
):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(
        mechanism,
        hformat,
        remove_1,
        remove_pow,
        min_op_count,
        recursive_op_count,
        store_in_jacobian,
        round_decimals,
    )
    conv.writer()
    conv.formatter()


def convert_lst(lst, hformat):
    """Convert mechanisms from a file containing a list of directories."""
    lpath = pathlib.Path(lst)
    with open(lst, "r") as f:
        for line in f:
            if not line.startswith("#"):
                mechname = lpath.parents[0] / line.strip()
                print(f"""Converting file {mechname}""")
                convert(mechname, hformat)


def main():
    """Convert cantera mechanisms to C++ files."""
    parser = argparse.ArgumentParser(description="Mechanism converter")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument(
        "-l", "--lst", help="Mechanism directory file list", type=str
    )

    parser.add_argument(
        "--hformat",
        help="sytle format for .H file output",
        type=str,
        choices=["readable", "gpu"],
        default="gpu",
        required=False,
    )

    parser.add_argument(
        "-r1",
        "--remove_1",
        action="store_true",
        help="Remove factor 1.0 in printed expressions",
    )

    parser.add_argument(
        "-rp",
        "--remove_pow",
        action="store_true",
        help="Remove pow(...,n) if n=2,3 in printed expressions",
    )

    parser.add_argument(
        "-moc",
        "--min_op_count",
        type=int,
        metavar="",
        required=False,
        help="Min number of operation count per expression",
        default=0,
    )

    parser.add_argument(
        "-roc",
        "--recursive_op_count",
        action="store_true",
        help="Recursive elimination of expression",
    )

    parser.add_argument(
        "-sj",
        "--store_in_jacobian",
        action="store_true",
        help="Store temporary arrays in Jacobian array",
    )

    parser.add_argument(
        "-rd",
        "--round_decimals",
        action="store_true",
        help="Round decimal numbers when possible",
    )

    args = parser.parse_args()

    if args.fname:
        convert(
            args.fname,
            args.hformat,
            args.remove_1,
            args.remove_pow,
            args.min_op_count,
            args.recursive_op_count,
            args.store_in_jacobian,
            args.round_decimals,
        )
    elif args.lst:
        convert_lst(
            args.lst,
            args.hformat,
            args.remove_1,
            args.remove_pow,
            args.min_op_count,
            args.recursive_op_count,
            args.store_in_jacobian,
            args.round_decimals,
        )


if __name__ == "__main__":
    main()

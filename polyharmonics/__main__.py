# type: ignore[attr-defined]
from enum import Enum
from random import choice
from typing import Optional

import typer
from rich.console import Console
from sympy import latex, pretty

from polyharmonics import legendre, version


class Color(str, Enum):
    white = "white"
    red1 = "red"
    cyan1 = "cyan"
    magenta1 = "magenta"
    yellow1 = "yellow"
    green1 = "green"


app = typer.Typer(
    name="polyharmonics",
    help="Ortogonal Polynomials in the unit sphere.",
    add_completion=False,
)
console = Console()


def version_callback(print_version: bool) -> None:
    """Print the version of the package."""
    if print_version:
        console.print(f"[yellow]polyharmonics[/] version: [bold blue]{version}[/]")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "-v",
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the polyharmonics package.",
    ),
):
    pass


@app.command(name="legendre")
def legendre_command(
    n: str = typer.Option(
        ...,
        help="""The degree of the polynomial(s).
        An integer or a comma-separated list of integers.""",
    ),
    print_latex: bool = typer.Option(
        False,
        "-l",
        "--latex",
        case_sensitive=False,
        help="Print the polynomial(s) in LaTeX format.",
    ),
    color: Optional[Color] = typer.Option(
        None,
        "-c",
        "--color",
        case_sensitive=False,
        help="Color for print. White if not specified.",
    ),
) -> None:
    """Calculate and print the Legendre polynomial(s)."""
    if color is None:
        color = Color.white

    # Convert the input to a list of integers
    try:
        n_values = [int(value) for value in n.split(",")]
        if any(i < 0 for i in n_values):
            raise typer.BadParameter("All integers must be greater or equal to 0")
    except ValueError:
        raise typer.BadParameter(
            "n must be an integer or a comma-separated list of integers"
        )

    # Calculate the Legendre polynomial(s)
    result = legendre(n_values)

    for n, pol in zip(n_values, result):
        console.print(f"[bold {color}]Legendre polynomial of degree {n}:")
        if print_latex:
            console.print(f"[bold {color}]{latex(pol)}[/]")
        else:
            console.print(f"[bold {color}]{pretty(pol)}[/]")


@app.command(name="assoc-legendre")
def assoc_legendre(
    color: Optional[Color] = typer.Option(
        None,
        "-c",
        "--color",
        case_sensitive=False,
        help="Color for print. If not specified then choice will be random.",
    ),
) -> None:
    """TODO"""
    if color is None:
        color = choice(list(Color))

    test: str = "Hi!"
    console.print(f"[bold {color}]{test}[/]")


if __name__ == "__main__":
    app()

#!/usr/bin/env python3

import argparse
import os
import sys
from typing import Optional, Tuple

import numpy as np


X_AXIS_CHOICES = ("kperp", "kpar", "kmag")
SCALE_CHOICES = ("linear", "log")
X_LABELS = {
    "kperp": r"$k_{\perp} d_R$",
    "kpar": r"$k_{\parallel} d_R$",
    "kmag": r"$|k| d_R$",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a two-panel quick-look dispersion plot from a four-column "
            "ALPS scan file."
        )
    )
    parser.add_argument("input_file", help="Path to the ALPS scan file to plot.")
    parser.add_argument(
        "--x-axis",
        choices=X_AXIS_CHOICES,
        default="kperp",
        help="Quantity to use on the x-axis.",
    )
    parser.add_argument(
        "--x-scale",
        choices=SCALE_CHOICES,
        default="linear",
        help="Scale to use for the shared x-axis.",
    )
    parser.add_argument(
        "--y-scale",
        choices=SCALE_CHOICES,
        default="linear",
        help="Scale to use for both y-axes.",
    )
    parser.add_argument(
        "--output",
        help="Output PNG path. Defaults to a name derived from the input file.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively after saving it.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output DPI for the saved figure.",
    )
    parser.add_argument(
        "--title",
        help="Optional figure title. Defaults to the input file basename.",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"Error: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_scan_file(path: str) -> np.ndarray:
    if not os.path.isfile(path):
        fail(f"input file does not exist: {path}")

    try:
        data = np.loadtxt(path)
    except Exception as exc:
        fail(f"could not read '{path}': {exc}")

    if data.ndim == 1:
        data = data[np.newaxis, :]

    if data.ndim != 2 or data.shape[1] < 4:
        fail(
            f"expected at least 4 numeric columns in '{path}', "
            f"found shape {data.shape}"
        )

    return data


def get_x_values(data: np.ndarray, x_axis: str) -> np.ndarray:
    kperp = data[:, 0].astype(float)
    kpar = data[:, 1].astype(float)

    if x_axis == "kperp":
        return kperp
    if x_axis == "kpar":
        return kpar
    if x_axis == "kmag":
        return np.sqrt(kperp**2 + kpar**2)
    fail(f"unsupported x-axis choice: {x_axis}")
    raise AssertionError("unreachable")


def validate_log_x(x_values: np.ndarray, x_label: str) -> None:
    non_positive = np.count_nonzero(x_values <= 0.0)
    if non_positive:
        fail(
            f"log x-axis requested for {x_label}, but found {non_positive} "
            "non-positive x values"
        )


def validate_log_y(values: np.ndarray, quantity_name: str) -> None:
    zero_count = np.count_nonzero(values == 0.0)
    if zero_count:
        fail(
            f"log y-axis requested for {quantity_name}, but found {zero_count} "
            "zero values that cannot be shown on a log axis"
        )


def sort_series(x_values: np.ndarray, y_values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(x_values)
    return x_values[order], y_values[order]


def plot_signed_log_line(ax, x_values: np.ndarray, y_values: np.ndarray) -> None:
    positive = np.where(y_values > 0.0, y_values, np.nan)
    negative = np.where(y_values < 0.0, np.abs(y_values), np.nan)

    if np.isfinite(positive).any():
        ax.plot(x_values, positive, linestyle="-", color="tab:blue", linewidth=1.8)
    if np.isfinite(negative).any():
        ax.plot(x_values, negative, linestyle="--", color="tab:blue", linewidth=1.8)


def plot_linear_line(ax, x_values: np.ndarray, y_values: np.ndarray) -> None:
    ax.plot(x_values, y_values, linestyle="-", color="tab:blue", linewidth=1.8)


def default_output_path(input_file: str, x_axis: str, x_scale: str, y_scale: str) -> str:
    basename = os.path.basename(input_file)
    filename = f"{basename}_{x_axis}_{x_scale}_{y_scale}.png"
    return os.path.join(os.getcwd(), filename)


def configure_axes(
    axes,
    x_label: str,
    x_scale: str,
    y_scale: str,
) -> None:
    for ax in axes:
        ax.set_xscale(x_scale)
        ax.set_yscale(y_scale)
        ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.5)

    axes[0].set_ylabel(r"$\omega_r/\Omega_R$")
    axes[1].set_ylabel(r"$\gamma/\Omega_R$")
    axes[1].set_xlabel(x_label)


def prepare_plot_data(
    data: np.ndarray,
    x_axis: str,
    x_scale: str,
    y_scale: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    x_values = get_x_values(data, x_axis)
    omega_r = data[:, 2].astype(float)
    gamma = data[:, 3].astype(float)
    x_label = X_LABELS[x_axis]

    if x_scale == "log":
        validate_log_x(x_values, x_label)

    if y_scale == "log":
        validate_log_y(omega_r, "omega_r")
        validate_log_y(gamma, "gamma")

    x_omega, omega_r = sort_series(x_values, omega_r)
    x_gamma, gamma = sort_series(x_values, gamma)
    return x_omega, omega_r, gamma, x_label


def draw_dispersion_axes(
    axes,
    data: np.ndarray,
    x_axis: str,
    x_scale: str,
    y_scale: str,
) -> str:
    x_values, omega_r, gamma, x_label = prepare_plot_data(data, x_axis, x_scale, y_scale)

    if y_scale == "log":
        plot_signed_log_line(axes[0], x_values, omega_r)
        plot_signed_log_line(axes[1], x_values, gamma)
    else:
        plot_linear_line(axes[0], x_values, omega_r)
        plot_linear_line(axes[1], x_values, gamma)

    configure_axes(list(axes), x_label, x_scale, y_scale)
    return x_label


def build_plot(
    data: np.ndarray,
    x_axis: str,
    x_scale: str,
    y_scale: str,
    title: Optional[str],
):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(9.0, 6.5),
        constrained_layout=True,
    )

    draw_dispersion_axes(axes, data, x_axis, x_scale, y_scale)

    if title:
        fig.suptitle(title)

    return fig


def main() -> int:
    args = parse_args()
    data = load_scan_file(args.input_file)

    title = args.title if args.title is not None else os.path.basename(args.input_file)
    output_path = args.output or default_output_path(
        args.input_file, args.x_axis, args.x_scale, args.y_scale
    )

    fig = build_plot(data, args.x_axis, args.x_scale, args.y_scale, title)
    import matplotlib.pyplot as plt

    try:
        fig.savefig(output_path, dpi=args.dpi)
    except Exception as exc:
        plt.close(fig)
        fail(f"could not save figure to '{output_path}': {exc}")

    print(f"Saved plot to {output_path}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

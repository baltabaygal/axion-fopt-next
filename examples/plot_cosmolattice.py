"""Plot a small CosmoLattice adapter run (requires matplotlib).

Example: python examples/plot_cosmolattice.py --run outputs/one_bubble_evolution
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_run(run, background, output):
    metadata = json.loads(background.read_text())
    n, length = metadata["n"], metadata["box_size"]
    history = np.genfromtxt(run / "history.csv", delimiter=",", names=True)
    final = np.genfromtxt(run / "final.csv", delimiter=",", names=True)
    if final.size != n**3:
        raise ValueError("Snapshot and background grid sizes disagree")
    np.testing.assert_array_equal(np.column_stack([final["x"],final["y"],final["z"]]),
                                 np.indices((n,n,n)).reshape(3,-1).T)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size":11, "axes.spines.top":False,
                        "axes.spines.right":False, "figure.dpi":160,
                        "savefig.dpi":200})
    mean = history["theta_mean"]
    variance = history["theta2_mean"] - mean**2
    # Moments from homogeneous fields can differ at roundoff level.
    variance[np.abs(variance) < 1e-13] = 0.
    if np.any(variance < -1e-13):
        raise ValueError("Negative field variance beyond roundoff")
    sigma = np.sqrt(np.maximum(variance, 0.))
    fig, axes = plt.subplots(2,1,figsize=(8,6),sharex=True,layout="constrained")
    axes[0].fill_between(history["tau"],mean-sigma,mean+sigma,
                         color="#407caa",alpha=.22,label=r"Spatial mean $\pm\,\sigma_\theta$")
    axes[0].plot(history["tau"],mean,color="#24587b",lw=2,label=r"$\langle\theta\rangle$")
    axes[0].axhline(0,color=".5",lw=.7)
    axes[0].set_ylabel(r"Axion angle $\theta=\phi/f$")
    axes[0].set_title(f"One-bubble axion evolution · {n}³ periodic lattice",loc="left",pad=12)
    axes[0].legend(frameon=False,fontsize=10)
    axes[1].plot(history["tau"],sigma,color="#c57238",lw=2)
    axes[1].set_ylabel(r"Spatial spread $\sigma_\theta$")
    axes[1].set_xlabel(r"Conformal time $\tau$ (simulation units)")
    for ax in axes:
        ax.grid(alpha=.18)
        ax.set_xlim(history["tau"][0],history["tau"][-1])
    fig.savefig(output / "field_evolution.png")
    fig.savefig(output / "field_evolution.pdf")
    plt.close(fig)

    birth_centre = metadata["centres"][0]
    z_index = int(np.floor(birth_centre[2] / length * n)) % n
    z = (z_index+.5)*length/n
    arrival = final["arrival"].reshape(n,n,n)[:,:,z_index]
    theta = final["theta"].reshape(n,n,n)[:,:,z_index]
    fig, axes = plt.subplots(1,2,figsize=(10,4.5),layout="constrained")
    extent = (0,length,0,length)
    first = axes[0].imshow(np.ma.masked_invalid(arrival.T),origin="lower",extent=extent,
                           cmap="viridis",interpolation="nearest")
    fig.colorbar(first,ax=axes[0],label=r"First-arrival time $\tau_{\rm wall}$",shrink=.85)
    axes[0].set_title("Prescribed bubble geometry",loc="left")
    axes[0].plot(birth_centre[0],birth_centre[1],marker="+",color="white",ms=9,mew=1.5)
    bound = float(np.max(np.abs(theta)))
    second = axes[1].imshow(theta.T,origin="lower",extent=extent,cmap="RdBu_r",
                            vmin=-bound,vmax=bound,interpolation="nearest")
    fig.colorbar(second,ax=axes[1],label=r"Final axion angle $\theta$",shrink=.85)
    axes[1].set_title(f"Evolved field at τ = {history['tau'][-1]:g}",loc="left")
    for ax in axes:
        ax.set_xlabel("Comoving x (simulation units)")
        ax.set_ylabel("Comoving y (simulation units)")
    fig.suptitle(f"One-bubble spatial slice · z = {z:g} · periodic box L = {length:g}",fontsize=12)
    fig.savefig(output / "bubble_slice.png")
    fig.savefig(output / "bubble_slice.pdf")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--background",type=Path,default=Path("outputs/one_bubble.json"))
    parser.add_argument("--output",type=Path,default=Path("outputs/plots"))
    args=parser.parse_args()
    plot_run(args.run,args.background,args.output)

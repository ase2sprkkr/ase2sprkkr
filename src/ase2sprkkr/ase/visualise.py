from typing import Tuple, Union
from ase.io import write
from io import BytesIO
import sys
import subprocess


def view(
    atoms,
    repeat: Union[int, Tuple[int, int, int]] = None,
    scale_radii: float = 0.5,
    rotations: str = "",
    bonds=False,
    graph: str = "",
    terminal=False,
    blocking: bool = False,
    close_after: float | None = None,
):
    """Visualise atoms using ase viewer

    Parameters
    ----------
    repeat
      Repeat the structure along axes [x,y,z]. If int is given, use the same number for all axis.

    scale_radii
      Relative size of atomic sites

    rotations
      Examples: "x90", "y90,z90".

    bonds
      Show bonds between atoms.

    graph

    terminal

    blocking
      If False, launch the viewer in a subprocess and return immediately.

    close_after
      Seconds after which the Tk viewer exits its event loop. When set,
      the viewer runs in the current process, even if blocking is False.
    """
    from ase.gui.images import Images

    if close_after is not None and close_after < 0:
        raise ValueError("close_after must not be negative")

    if repeat is not None:
        if isinstance(repeat, int):
            repeat = [repeat] * 3

    if not blocking and close_after is None:
        buf = BytesIO()
        write(buf, atoms, format="traj")

        args = [sys.executable, "-m", "ase", "gui", "-"]
        if repeat:
            args.append("--repeat={},{},{}".format(*repeat))
        if bonds:
            args.append("-b")
        if scale_radii:
            args.append("--scale={}".format(scale_radii))
        if rotations:
            args.append("--rotations={}".format(rotations))
        if graph:
            args.append("--graph={}".format(graph))
        if terminal:
            args.append("--terminal={}".format(terminal))

        proc = subprocess.Popen(args,
                                stdin=subprocess.PIPE,
                                start_new_session=True,
                                )
        proc.stdin.write(buf.getvalue())
        proc.stdin.close()
        return proc

    images = Images()
    images.initialize([atoms])

    if repeat:
        images.repeat_images(repeat)
    if scale_radii:
        images.scale_radii(scale_radii)

    if terminal:
        if args.graph is not None:
            data = images.graph(args.graph)
            for line in data.T:
                for x in line:
                    print(x, end=" ")
                print()
    else:
        import os
        from ase.gui.gui import GUI

        backend = os.environ.get("MPLBACKEND", "")
        if backend == "module://ipykernel.pylab.backend_inline":
            # Jupyter should not steal our windows
            del os.environ["MPLBACKEND"]

        gui = GUI(images, rotations, bonds, graph)
        if close_after is not None:
            gui.window.win.after(round(close_after * 1000), gui.window.win.quit)
        gui.run()

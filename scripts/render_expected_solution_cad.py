#!/usr/bin/env python3
"""Build and render the expected two-gate proof architecture as parametric CAD.

Outputs a STEP assembly and a deterministic headless PNG.
The model is explanatory rather than a claim about a literal geometric limit.
"""
from __future__ import annotations

import math
from pathlib import Path

import cadquery as cq
import vtk


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "evidence" / "share"
RADIUS = 40.0
SURFACE = RADIUS + 0.8


def box(x0: float, x1: float, width: float = 90.0) -> cq.Shape:
    return cq.Solid.makeBox(
        x1 - x0,
        width,
        width,
        (x0, -width / 2, -width / 2),
    )


def rod(first: tuple[float, float, float], second: tuple[float, float, float], radius: float) -> cq.Shape:
    vector = cq.Vector(*(b - a for a, b in zip(first, second)))
    return cq.Solid.makeCylinder(radius, vector.Length, first, vector.normalized())


def conic_nodes() -> tuple[tuple[float, float, float], ...]:
    """A deterministic spherical projection of t -> (t^2,t,1), t in F_17."""
    nodes = []
    for t in range(16):
        z = SURFACE * 0.72 * (2 * t / 15 - 1)
        phase = (((t * t + 3 * t + 1) % 17) / 16 - 0.5) * 1.95
        radial = math.sqrt(SURFACE * SURFACE - z * z)
        nodes.append((radial * math.cos(phase), radial * math.sin(phase), z))
    return tuple(nodes)


def build_parts() -> list[tuple[str, cq.Shape, tuple[float, float, float], float]]:
    outer = cq.Solid.makeSphere(RADIUS, angleDegrees1=-90, angleDegrees2=90)
    inner = cq.Solid.makeSphere(RADIUS - 1.1, angleDegrees1=-90, angleDegrees2=90)
    shell = outer.cut(inner)
    notch = cq.Solid.makeBox(7, 90, 13, (-3.5, -45, 28))
    left = shell.intersect(box(-43, 0)).cut(notch)
    right = shell.intersect(box(0, 43)).cut(notch)

    parts: list[tuple[str, cq.Shape, tuple[float, float, float], float]] = [
        ("spectral shell", left, (0.05, 0.38, 0.62), 0.34),
        ("incidence shell", right, (0.66, 0.34, 0.035), 0.34),
    ]

    left_clip = box(-44, -0.15, 94)
    for index, z in enumerate((-29, -20, -11, 0, 11, 20, 29)):
        major = math.sqrt((SURFACE + 0.25) ** 2 - z * z)
        ring = cq.Solid.makeTorus(major, 0.34, (0, 0, z)).intersect(left_clip)
        parts.append((f"spectral latitude {index}", ring, (0.05, 0.72, 1.0), 1.0))
    for index, angle in enumerate((-52, -24, 24, 52)):
        radians = math.radians(angle)
        axis = (0, math.cos(radians), math.sin(radians))
        ring = cq.Solid.makeTorus(SURFACE + 0.25, 0.30, dir=axis).intersect(left_clip)
        parts.append((f"spectral meridian {index}", ring, (0.10, 0.56, 0.96), 1.0))

    right_clip = box(0.15, 44, 94)
    for index, z in enumerate((-28, -14, 0, 14, 28)):
        major = math.sqrt((SURFACE + 0.18) ** 2 - z * z)
        ring = cq.Solid.makeTorus(major, 0.20, (0, 0, z)).intersect(right_clip)
        parts.append((f"incidence latitude {index}", ring, (0.66, 0.31, 0.035), 0.62))
    for index, angle in enumerate((-42, 0, 42)):
        radians = math.radians(angle)
        axis = (0, math.cos(radians), math.sin(radians))
        ring = cq.Solid.makeTorus(SURFACE + 0.18, 0.18, dir=axis).intersect(right_clip)
        parts.append((f"incidence meridian {index}", ring, (0.55, 0.25, 0.025), 0.58))

    seam = cq.Solid.makeTorus(RADIUS + 0.2, 0.72, dir=(1, 0, 0)).cut(notch)
    parts.append(("implication seam", seam, (0.92, 0.95, 1.0), 1.0))

    nodes = conic_nodes()
    for index, point in enumerate(nodes):
        marker = cq.Solid.makeSphere(0.95, point, angleDegrees1=-90, angleDegrees2=90)
        parts.append((f"finite point {index}", marker, (1.0, 0.62, 0.08), 1.0))
    # Polygonal conic trace and four distinguished secants.
    for index in range(len(nodes) - 1):
        parts.append(
            (f"conic trace {index}", rod(nodes[index], nodes[index + 1], 0.23), (0.95, 0.40, 0.04), 1.0)
        )
    for index, (first, second) in enumerate(((0, 8), (2, 11), (4, 14), (6, 15))):
        parts.append(
            (f"secant {index}", rod(nodes[first], nodes[second], 0.34), (1.0, 0.78, 0.18), 1.0)
        )
    return parts


def vtk_actor(shape: cq.Shape, color: tuple[float, float, float], opacity: float):
    vertices, triangles = shape.tessellate(0.28, 0.12)
    points = vtk.vtkPoints()
    for vertex in vertices:
        points.InsertNextPoint(vertex.x, vertex.y, vertex.z)
    cells = vtk.vtkCellArray()
    for triangle in triangles:
        cell = vtk.vtkTriangle()
        for index, point_index in enumerate(triangle):
            cell.GetPointIds().SetId(index, point_index)
        cells.InsertNextCell(cell)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.ConsistencyOn()
    normals.SplittingOff()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().SetAmbient(0.28)
    actor.GetProperty().SetDiffuse(0.72)
    actor.GetProperty().SetSpecular(0.45)
    actor.GetProperty().SetSpecularPower(32)
    return actor


def render(parts, target: Path) -> None:
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.006, 0.009, 0.016)
    renderer.SetUseDepthPeeling(False)
    for _name, shape, color, opacity in parts:
        renderer.AddActor(vtk_actor(shape, color, opacity))

    key = vtk.vtkLight()
    key.SetPosition(-80, -120, 100)
    key.SetFocalPoint(0, 0, 0)
    key.SetIntensity(0.9)
    renderer.AddLight(key)
    fill = vtk.vtkLight()
    fill.SetPosition(90, -70, 20)
    fill.SetFocalPoint(0, 0, 0)
    fill.SetIntensity(0.65)
    renderer.AddLight(fill)

    camera = renderer.GetActiveCamera()
    camera.SetPosition(0, -155, 13)
    camera.SetFocalPoint(0, 0, 0)
    camera.SetViewUp(0, 0, 1)
    camera.ParallelProjectionOn()
    camera.SetParallelScale(50)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(True)
    window.SetAlphaBitPlanes(True)
    window.SetMultiSamples(8)
    window.SetSize(1536, 900)
    window.AddRenderer(renderer)
    window.Render()
    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetScale(1)
    capture.ReadFrontBufferOff()
    capture.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(target))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    parts = build_parts()
    compound = cq.Compound.makeCompound([shape for _name, shape, _color, _opacity in parts])
    cq.exporters.export(compound, str(OUTPUT / "expected-solution-structure.step"))
    render(parts, OUTPUT / "expected-solution-structure.png")
    print(f"wrote CAD assembly and render to {OUTPUT}")


if __name__ == "__main__":
    main()

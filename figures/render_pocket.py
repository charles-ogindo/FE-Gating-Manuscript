#!/usr/bin/env python3
"""Regenerate Figure 1 (paclitaxel in the beta-tubulin taxane pocket) with PyMOL.

Usage (from anywhere):
    pip install pymol-open-source
    python figures/render_pocket.py

Reads  structures/complex.pdb  and writes  figures/paclitaxel_pocket_render.png.

The pocket residues are selected by their position in complex.pdb (which renumbers
the chain) and labelled with the manuscript's beta-tubulin author numbering, matched
by residue identity and proximity to the ligand. Hydrogen bonds (Thr276, Gly370) are
drawn as measured dashed lines; the salt bridge (Asp26) and pi-cation (His229) appear
as labelled residues. Geometry and interaction identities are taken from the pose
contact artifact contacts_00.json (schema 1.3.0).
"""
import os
import pymol

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
COMPLEX = os.path.join(ROOT, "structures", "complex.pdb")
OUT = os.path.join(HERE, "paclitaxel_pocket_render.png")

pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd, util

cmd.load(COMPLEX, "cx")
cmd.hide("everything")
cmd.select("lig", "resn LIG")

# file residue number -> manuscript (author) label, matched by identity + proximity
m = {22: "Val23", 25: "Asp26", 226: "His229", 271: "Pro274",
     273: "Thr276", 359: "Gly370", 360: "Leu371"}
psel = "chain B and resi " + "+".join(str(k) for k in m)
cmd.select("pocket", psel)

# beta-tubulin cartoon, soft and transparent for context
cmd.show("cartoon", "chain B and polymer")
cmd.set("cartoon_transparency", 0.62, "chain B")
cmd.color("gray80", "chain B and polymer")

# ligand
cmd.show("sticks", "lig")
cmd.set("stick_radius", 0.19, "lig")
cmd.color("yelloworange", "lig and elem C")
util.cnc("lig")

# pocket residues
cmd.show("sticks", "pocket")
cmd.set("stick_radius", 0.15, "pocket")
cmd.color("palecyan", "pocket and elem C")
util.cnc("pocket")

cmd.hide("everything", "hydro")  # drop hydrogens for clarity

# hydrogen bonds (polar contacts) ligand <-> pocket
cmd.distance("hb", "lig", "pocket and (elem N+O)", mode=2)
cmd.color("grey20", "hb")
cmd.set("dash_width", 3.2)
cmd.set("dash_gap", 0.32)
cmd.hide("labels", "hb")

# author-number labels
for r, lab in m.items():
    cmd.label(f"chain B and resi {r} and name CA", f'"{lab}"')
cmd.set("label_size", 16)
cmd.set("label_color", "black")
cmd.set("label_font_id", 7)

# quality and view
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 1)
cmd.set("antialias", 2)
cmd.set("ambient", 0.4)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("specular", 0.2)
cmd.orient("lig or pocket")
cmd.zoom("lig or pocket", 3.2)
cmd.turn("y", 12)
cmd.turn("x", -5)
cmd.ray(1500, 1150)
cmd.png(OUT, dpi=200)
print("wrote", OUT)

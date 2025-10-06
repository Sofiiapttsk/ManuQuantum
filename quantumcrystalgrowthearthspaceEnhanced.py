import sys
import math
import random
import numpy as np

from PyQt5.QtCore import Qt, QTimer, QPointF
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QSlider, QLabel, QComboBox, QPushButton,
                             QGroupBox, QRadioButton, QTabWidget, QCheckBox, QFrame)
from PyQt5.QtGui import QColor, QPainter, QBrush, QPen, QFont
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from PyQt5.QtOpenGL import QGLWidget

# DEBUG  

DEBUG_MODE = False  # set True to get console prints

def dprint(*a):
    if DEBUG_MODE:
        print("[DEBUG]", *a)

# PHYSICAL MODEL PARAMETERS

RANDOM_SEED = 1234
random.seed(RANDOM_SEED)

DOPANT_N = 64
DOPANT_RADIUS = 1.0
DOPANT_HEIGHT = 1.2
BASE_DIFFUSIVITY_SPACE = 2.5e-3
BASE_DIFFUSIVITY_EARTH = 1.2e-3
ADVECT_SCALE = 0.6
MAX_DT = 0.05
INTERFACE_THICKNESS = 0.08
INCORPORATION_RATE = 0.15
UNIFORMITY_SMOOTH = 0.1

ACOUSTIC_WAVE_NUMBER = 6.0
ACOUSTIC_DAMP = 0.9
ACOUSTIC_SURFACE_MOD = 0.03
ACOUSTIC_STABILIZATION = 0.25
ACOUSTIC_UNIFORMITY_BOOST = 0.08

QUALITY_WEIGHT_DOPANT_UNIFORM = 0.35
QUALITY_WEIGHT_TEMP = 0.30
QUALITY_WEIGHT_ACOUSTIC = 0.15
QUALITY_WEIGHT_DEFECT = 0.20


# MATERIALS 

QUANTUM_MATERIALS = {
    "Topological Insulator (Bi₂Se₃)": {
        "melting_point": 705,
        "crystal_structure": "Rhombohedral",
        "growth_rate": 0.8,
        "color": (0.2, 0.4, 0.8),
        "defect_earth": 0.15,
        "defect_space": 0.02,
        "thermal_conductivity": 2.8,
        "dopant_diffusivity_earth": 0.4,
        "dopant_diffusivity_space": 0.9,
        "typical_defects": ["Vacancies", "Dislocations", "Stacking faults"],
        "space_advantages": ["No container wall interactions", "Uniform temperature gradient"]
    },
    "High-Tc Superconductor (YBCO)": {
        "melting_point": 1200,
        "crystal_structure": "Orthorhombic",
        "growth_rate": 0.5,
        "color": (0.1, 0.7, 0.2),
        "defect_earth": 0.2,
        "defect_space": 0.03,
        "thermal_conductivity": 5.7,
        "dopant_diffusivity_earth": 0.3,
        "dopant_diffusivity_space": 0.85,
        "typical_defects": ["Oxygen vacancies", "Grain boundaries", "Twin boundaries"],
        "space_advantages": ["Homogeneous oxygen distribution", "Reduced phase separation"]
    },
    "Quantum Spin Liquid (α-RuCl₃)": {
        "melting_point": 850,
        "crystal_structure": "Honeycomb Lattice",
        "growth_rate": 0.3,
        "color": (0.8, 0.2, 0.5),
        "defect_earth": 0.25,
        "defect_space": 0.04,
        "thermal_conductivity": 4.2,
        "dopant_diffusivity_earth": 0.35,
        "dopant_diffusivity_space": 0.95,
        "typical_defects": ["Honeycomb lattice distortions", "Stacking faults", "Impurity phases"],
        "space_advantages": ["Perfect honeycomb arrangement", "Reduced strain at layer interfaces"]
    },
}

# Safe fallback material if a key mismatch occurs 
FALLBACK_MATERIAL_KEY = list(QUANTUM_MATERIALS.keys())[0]


# Utility

def clamp(v, lo, hi):
    return lo if v < lo else hi if v > hi else v

def compute_normal(v1, v2, v3):
    ux, uy, uz = v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]
    vx, vy, vz = v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2]
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    length = math.sqrt(nx*nx + ny*ny + nz*nz)
    if length > 1e-9:
        return nx/length, ny/length, nz/length
    return 0.0, 0.0, 1.0

def get_material_data(name):
    """Robust accessor (PATCH)"""
    if name in QUANTUM_MATERIALS:
        return QUANTUM_MATERIALS[name]
    dprint("Material key missing:", name, "falling back to", FALLBACK_MATERIAL_KEY)
    return QUANTUM_MATERIALS[FALLBACK_MATERIAL_KEY]


# 2D Convection Visualizer

class ConvectionVisualizer(QWidget):
    def __init__(self, simulation_ref, parent=None):
        super(ConvectionVisualizer, self).__init__(parent)
        self.simulation = simulation_ref
        self.setMinimumSize(320, 320)
        self.particles = []
        self.dopants = []
        self.show_convection = True
        self.show_dopants = True
        self.show_temperature = True
        self.reset_particles()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_particles)
        self.timer.start(50)

    def reset_particles(self):
        self.particles.clear()
        self.dopants.clear()
        for _ in range(80):
            self.particles.append({
                'x': random.uniform(0.05, 0.95),
                'y': random.uniform(0.05, 0.95),
                'size': random.uniform(2, 4)
            })
        for _ in range(40):
            self.dopants.append({
                'x': random.uniform(0.1, 0.9),
                'y': random.uniform(0.1, 0.9),
                'size': random.uniform(3, 6)
            })

    def _velocity_field(self, x, y):
        if not self.simulation or self.simulation.gravity <= 0.5:
            return 0.0, 0.0
        xm = (x - 0.5)
        ym = (y - 0.5)
        dpsi_dx = math.pi * math.cos(math.pi * xm) * math.sin(math.pi * ym)
        dpsi_dy = math.pi * math.sin(math.pi * xm) * math.cos(math.pi * ym)
        vr = dpsi_dy
        vt = -dpsi_dx
        r2 = xm*xm + ym*ym
        damp = (1.0 - min(1.0, r2 * 4.0))
        acoustic_factor = 1.0 - ACOUSTIC_STABILIZATION * self.simulation.acoustic_intensity
        scale = 0.0025 * damp * acoustic_factor
        return vr * scale, vt * scale

    def update_particles(self):
        # guard simulation existence + try/except to catch key errors cleanly
        if not self.isVisible():
            return
        if not hasattr(self.simulation, "material"):
            return
        try:
            g = self.simulation.gravity
        except Exception as e:
            dprint("Accessing simulation failed:", e)
            return

        for p in self.particles:
            try:
                if g > 0.5:
                    vx, vy = self._velocity_field(p['x'], p['y'])
                    p['x'] += vx + random.uniform(-0.0002, 0.0002)
                    p['y'] += vy + random.uniform(-0.0002, 0.0002)
                else:
                    p['x'] += random.uniform(-0.0004, 0.0004)
                    p['y'] += random.uniform(-0.0004, 0.0004)
            except Exception as e:
                dprint("Particle update error:", e)
            p['x'] = (p['x'] + 1.0) % 1.0
            p['y'] = (p['y'] + 1.0) % 1.0

        # Material lookup with fallback 
        mat = get_material_data(getattr(self.simulation, "material", FALLBACK_MATERIAL_KEY))

        for d in self.dopants:
            try:
                if g > 0.5:
                    vx, vy = self._velocity_field(d['x'], d['y'])
                    diff = 0.0004 * mat['dopant_diffusivity_earth']
                    d['x'] += vx + random.uniform(-diff, diff)
                    d['y'] += vy + random.uniform(-diff, diff)
                else:
                    diff = 0.0005 * mat['dopant_diffusivity_space']
                    d['x'] += random.uniform(-diff, diff)
                    d['y'] += random.uniform(-diff, diff)
                d['x'] = max(0.02, min(0.98, d['x']))
                d['y'] = max(0.02, min(0.98, d['y']))
            except Exception as e:
                dprint("Dopant update error:", e)
        self.update()

    def paintEvent(self, event):
        # Robustness:if simulation not ready - skip
        if not hasattr(self.simulation, "growth_progress"):
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        rect = self.rect()
        cx = rect.width() / 2
        cy = rect.height() / 2
        radius = min(cx, cy) * 0.85

        if self.show_temperature:
            for y in range(rect.height()):
                t = 1.0 - y / rect.height()
                if getattr(self.simulation, "gravity", 0.0) > 0.5:
                    red = int(255 * min(1.0, 0.5 + t * 0.9))
                    green = int(180 * min(1.0, 0.3 + t * 0.5))
                    blue = int(50 * t)
                else:
                    red = int(255 * min(1.0, 0.5 + t * 0.5))
                    green = int(150 * min(1.0, 0.3 + t * 0.5))
                    blue = int(50 * t)
                painter.setPen(QColor(red, green, blue))
                painter.drawLine(0, y, rect.width(), y)

        painter.setPen(QPen(QColor(200, 200, 200), 2))
        painter.setBrush(Qt.NoBrush)
        painter.drawEllipse(cx - radius, cy - radius, radius * 2, radius * 2)

        melt_radius = radius * (1.0 - self.simulation.growth_progress * 0.7)
        painter.setBrush(QBrush(QColor(255, 160, 20, 180)))
        painter.setPen(Qt.NoPen)
        painter.drawEllipse(cx - melt_radius, cy - melt_radius, melt_radius * 2, melt_radius * 2)

        if self.show_convection:
            painter.setPen(Qt.NoPen)
            for p in self.particles:
                x = cx + (p['x'] - 0.5) * radius * 2
                y = cy + (p['y'] - 0.5) * radius * 2
                col = QColor(210, 230, 255, 140) if self.simulation.gravity > 0.5 else QColor(200, 220, 255, 100)
                painter.setBrush(QBrush(col))
                painter.drawEllipse(QPointF(x, y), p['size'], p['size'])

        if self.show_dopants:
            painter.setPen(Qt.NoPen)
            for dpt in self.dopants:
                x = cx + (dpt['x'] - 0.5) * radius * 2
                y = cy + (dpt['y'] - 0.5) * radius * 2
                col = QColor(50, 200, 50, 190) if self.simulation.gravity > 0.5 else QColor(90, 230, 90, 190)
                painter.setBrush(QBrush(col))
                painter.drawEllipse(QPointF(x, y), dpt['size'], dpt['size'])

        if self.simulation.growth_progress > 0:
            crystal_radius = radius * self.simulation.growth_progress * 0.8
            mat = get_material_data(self.simulation.material)
            col = mat["color"]
            q = self.simulation.crystal_quality
            r = int(col[0] * 255 * (0.65 + 0.35 * q))
            g = int(col[1] * 255 * (0.65 + 0.35 * q))
            b = int(col[2] * 255 * (0.65 + 0.35 * q))
            painter.setBrush(QBrush(QColor(r, g, b, 210)))
            painter.setPen(QPen(QColor(r//2, g//2, b//2, 180), 2))
            painter.drawEllipse(QPointF(cx, cy), crystal_radius, crystal_radius)

        painter.setPen(QPen(QColor(255, 255, 255), 1))
        painter.setFont(QFont("Arial", 10, QFont.Bold))
        title = ("Earth: Structured Convection & Dopant Advection"
                 if self.simulation.gravity > 0.5 else
                 "Space: Diffusion-Dominated Transport")
        painter.drawText(10, 20, title)


# 3D Simulation Widget (unchanged core + small visibility guard PATCH)

class CrystalGrowthSimulation(QGLWidget):
    def __init__(self, parent=None):
        super(CrystalGrowthSimulation, self).__init__(parent)
        self.setMinimumSize(800, 600)
        self.gravity = 0.0
        self.growth_progress = 0.0
        self.melt_wobble = 0.0
        self.melt_scale = 1.0
        self.melt_color = (1.0, 0.5, 0.0)
        self.material = list(QUANTUM_MATERIALS.keys())[0]
        self.crystal_quality = 1.0
        self.acoustic_intensity = 0.5
        self.temperature = 800
        self.rotation_speed = 0.5
        self.view_angle = 0.0
        self.show_flow_lines = True
        self.show_dopants = True
        self.show_crystal_lattice = True
        self.show_defects = True
        self.auto_grow = False
        self.auto_grow_speed = 0.001
        self.animation_direction = 1

        self.acoustic_sources = [
            (-1.5, 0, 0), (1.5, 0, 0),
            (0, -1.5, 0), (0, 1.5, 0),
            (0, 0, -1.5), (0, 0, 1.5)
        ]
        self.crystal_vertices = []
        self.crystal_faces = []
        self.generate_crystal_structure()
        self.defects = []
        self.dopants = []
        self.generate_dopants()
        self.flow_lines = []
        self.generate_flow_lines()

        self.dopant_field = np.zeros((DOPANT_N, DOPANT_N), dtype=np.float32)
        self.dopant_field_next = np.zeros_like(self.dopant_field)
        self._init_dopant_field()
        self.total_incorporated_dopant = 0.0
        self.incorporated_positions = []
        self.displayed_uniformity = 1.0

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start(16)
        glutInit()

    def _init_dopant_field(self):
        for i in range(DOPANT_N):
            for j in range(DOPANT_N):
                r = i / (DOPANT_N - 1)
                self.dopant_field[i, j] = 1.0 - 0.5 * r

    def generate_crystal_structure(self):
        material_data = get_material_data(self.material)
        structure = material_data["crystal_structure"]
        self.crystal_vertices = []
        self.crystal_faces = []
        if structure == "Rhombohedral":
            r = 0.5
            h = 1.0
            for i in range(6):
                ang = 2 * math.pi * i / 6
                self.crystal_vertices.append((r * math.cos(ang), h / 2, r * math.sin(ang)))
                self.crystal_vertices.append((r * math.cos(ang), -h / 2, r * math.sin(ang)))
            self.crystal_faces.append([0, 2, 4, 6, 8, 10])
            self.crystal_faces.append([1, 3, 5, 7, 9, 11])
            for i in range(6):
                i1 = i * 2
                i2 = (i * 2 + 2) % 12
                i3 = (i * 2 + 3) % 12
                i4 = i * 2 + 1
                self.crystal_faces.append([i1, i2, i3, i4])
            if self.show_crystal_lattice:
                for layer in range(1, 3):
                    y = h / 2 - layer * h / 3
                    layer_verts = []
                    for i in range(6):
                        ang = 2 * math.pi * i / 6
                        scale = 0.9 - layer * 0.1
                        self.crystal_vertices.append(
                            (r * math.cos(ang) * scale, y, r * math.sin(ang) * scale)
                        )
                        layer_verts.append(len(self.crystal_vertices) - 1)
                    self.crystal_faces.append(layer_verts)
        elif structure == "Orthorhombic":
            self.crystal_vertices = [
                (0.5, 0.6, 0.4), (0.5, 0.6, -0.4),
                (-0.5, 0.6, -0.4), (-0.5, 0.6, 0.4),
                (0.5, -0.6, 0.4), (0.5, -0.6, -0.4),
                (-0.5, -0.6, -0.4), (-0.5, -0.6, 0.4),
            ]
            self.crystal_faces = [
                [0, 1, 2, 3],
                [4, 5, 6, 7],
                [0, 1, 5, 4],
                [2, 3, 7, 6],
                [0, 3, 7, 4],
                [1, 2, 6, 5],
            ]
            if self.show_crystal_lattice:
                for y in [-0.3, 0.0, 0.3]:
                    base = len(self.crystal_vertices)
                    scale = 0.8
                    self.crystal_vertices.extend([
                        (0.5 * scale, y, 0.4 * scale),
                        (0.5 * scale, y, -0.4 * scale),
                        (-0.5 * scale, y, -0.4 * scale),
                        (-0.5 * scale, y, 0.4 * scale)
                    ])
                    self.crystal_faces.append([base, base + 1, base + 2, base + 3])
        elif structure == "Honeycomb Lattice":
            layers = 5 if self.show_crystal_lattice else 3
            layer_spacing = 0.2
            r = 0.5
            for layer in range(layers):
                z = layer * layer_spacing - (layers - 1) * layer_spacing / 2
                base_idx = len(self.crystal_vertices)
                for i in range(6):
                    ang = 2 * math.pi * i / 6
                    self.crystal_vertices.append((r * math.cos(ang), z, r * math.sin(ang)))
                self.crystal_faces.append([base_idx + i for i in range(6)])
                if layer > 0:
                    prev = base_idx - 6
                    for i in range(6):
                        self.crystal_faces.append([
                            prev + i,
                            prev + (i + 1) % 6,
                            base_idx + (i + 1) % 6,
                            base_idx + i
                        ])
                if self.show_crystal_lattice and layer % 2 == 0:
                    center_base = len(self.crystal_vertices)
                    small_r = r * 0.4
                    for i in range(6):
                        ang = 2 * math.pi * i / 6
                        self.crystal_vertices.append((small_r * math.cos(ang), z, small_r * math.sin(ang)))
                    self.crystal_faces.append([center_base + i for i in range(6)])
                    for i in range(6):
                        self.crystal_faces.append([
                            base_idx + i,
                            base_idx + (i + 1) % 6,
                            center_base + (i + 1) % 6,
                            center_base + i
                        ])
        # Regenerate defects snapshot with new geometry context
        self.generate_defects()

    def generate_defects(self):
        self.defects = []
        mat = get_material_data(self.material)
        prob = mat["defect_space" if self.gravity < 0.5 else "defect_earth"]
        max_defects = int(100 * self.growth_progress * prob)
        for _ in range(max_defects):
            if self.gravity > 0.5:
                if random.random() < 0.7:
                    side = random.choice(['bottom', 'edge'])
                    if side == 'bottom':
                        x = random.uniform(-0.5, 0.5) * self.growth_progress
                        y = random.uniform(-0.6, -0.3) * self.growth_progress
                        z = random.uniform(-0.5, 0.5) * self.growth_progress
                    else:
                        ang = random.uniform(0, 2 * math.pi)
                        dist = random.uniform(0.8, 1.0) * 0.5
                        x = math.cos(ang) * dist * self.growth_progress
                        y = random.uniform(-0.6, 0.6) * self.growth_progress
                        z = math.sin(ang) * dist * self.growth_progress
                else:
                    x = random.uniform(-0.5, 0.5) * self.growth_progress
                    y = random.uniform(-0.6, 0.6) * self.growth_progress
                    z = random.uniform(-0.5, 0.5) * self.growth_progress
            else:
                x = random.uniform(-0.5, 0.5) * self.growth_progress
                y = random.uniform(-0.6, 0.6) * self.growth_progress
                z = random.uniform(-0.5, 0.5) * self.growth_progress
            size = random.uniform(0.01, 0.04)
            defect_type = random.choice(['point', 'line', 'cluster'])
            self.defects.append((x, y, z, size, defect_type))

    def generate_dopants(self):
        self.dopants = []
        count = int(200 * self.growth_progress)
        mat = get_material_data(self.material)
        for _ in range(count):
            if self.gravity > 0.5:
                if random.random() < 0.7:
                    y = random.choice([
                        random.uniform(-0.6, -0.4),
                        random.uniform(-0.2, 0.0),
                        random.uniform(0.2, 0.4)
                    ]) * self.growth_progress
                    ang = random.uniform(0, 2 * math.pi)
                    dist = random.uniform(0.1, 0.4)
                    x = math.cos(ang) * dist * self.growth_progress
                    z = math.sin(ang) * dist * self.growth_progress
                else:
                    x = random.uniform(-0.5, 0.5) * self.growth_progress
                    y = random.uniform(-0.6, 0.6) * self.growth_progress
                    z = random.uniform(-0.5, 0.5) * self.growth_progress
            else:
                x = random.uniform(-0.45, 0.45) * self.growth_progress
                y = random.uniform(-0.55, 0.55) * self.growth_progress
                z = random.uniform(-0.45, 0.45) * self.growth_progress
            size = random.uniform(0.01, 0.02)
            self.dopants.append((x, y, z, size))

    def generate_flow_lines(self):
        self.flow_lines = []
        if self.gravity < 0.5:
            return
        for _ in range(35):
            ang = random.uniform(0, 2 * math.pi)
            h = random.uniform(-0.8, 0.8)
            radius = 1.0
            x = math.cos(ang) * radius
            y = h
            z = math.sin(ang) * radius
            if math.sqrt(x * x + z * z) < 0.5:
                direction = (0, 1, 0)
            else:
                direction = (0, -1, 0)
            mag = math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2)
            if mag > 1e-9:
                direction = (direction[0]/mag, direction[1]/mag, direction[2]/mag)
            length = random.uniform(0.25, 0.55)
            color = (0.2, 0.4, 0.8, 0.45)
            self.flow_lines.append((x, y, z, direction, length, color))

    # Acoustic modulation reused from prior version (omitted for brevity)
    def _surface_acoustic_modulation(self, theta, phi, t):
        a = ACOUSTIC_SURFACE_MOD * self.acoustic_intensity
        y20 = 0.5 * (3*math.cos(phi)**2 - 1)
        y22 = math.sin(phi)**2 * math.cos(2*theta)
        temporal = math.sin(2*t)
        return a * (0.6*y20 + 0.4*y22) * temporal

    def _dopant_diffusivity(self):
        mat = get_material_data(self.material)
        if self.gravity < 0.5:
            return BASE_DIFFUSIVITY_SPACE * mat["dopant_diffusivity_space"]
        return BASE_DIFFUSIVITY_EARTH * mat["dopant_diffusivity_earth"]

    def _dopant_velocity_field(self, r, z):
        if self.gravity < 0.5:
            return 0.0, 0.0
        H = DOPANT_HEIGHT
        zn = (z + H/2)/H
        if zn < 0.0 or zn > 1.0 or r < 0.0 or r > 1.0:
            return 0.0, 0.0
        dpsi_dz = math.pi * math.cos(math.pi * zn) * r * (1 - r) / H
        dpsi_dr = math.sin(math.pi * zn) * (1 - 2*r)
        if r < 1e-4:
            v_r = 0.0
            v_z = -dpsi_dr / 0.5
        else:
            v_r = dpsi_dz / (r + 1e-6)
            v_z = -dpsi_dr / (r + 1e-6)
        acoustic_factor = 1.0 - ACOUSTIC_STABILIZATION * self.acoustic_intensity
        v_r *= ADVECT_SCALE * acoustic_factor
        v_z *= ADVECT_SCALE * acoustic_factor
        return v_r, v_z

    def _dopant_step(self, dt):
        D = self._dopant_diffusivity()
        f = self.dopant_field
        fn = self.dopant_field_next
        n = DOPANT_N
        dr = 1.0/(n-1)
        dz = DOPANT_HEIGHT/(n-1)
        for i in range(n):
            r = i*dr
            for j in range(n):
                z = j*dz - DOPANT_HEIGHT/2
                vr, vz = self._dopant_velocity_field(r, z)
                r_back = clamp(r - vr * dt, 0.0, 1.0)
                z_back = clamp(z - vz * dt, -DOPANT_HEIGHT/2, DOPANT_HEIGHT/2)
                ib = r_back / dr
                jb = (z_back + DOPANT_HEIGHT/2)/dz
                i0 = int(clamp(math.floor(ib), 0, n-2))
                j0 = int(clamp(math.floor(jb), 0, n-2))
                di = ib - i0
                dj = jb - j0
                adv_val = (f[i0,j0]*(1-di)*(1-dj) +
                           f[i0+1,j0]*di*(1-dj) +
                           f[i0,j0+1]*(1-di)*dj +
                           f[i0+1,j0+1]*di*dj)
                im = i-1 if i>0 else i
                ip = i+1 if i<n-1 else i
                jm = j-1 if j>0 else j
                jp = j+1 if j<n-1 else j
                lap = (f[ip,j] - 2*f[i,j] + f[im,j])/(dr*dr) + (f[i,jp] - 2*f[i,j] + f[i,jm])/(dz*dz)
                fn[i,j] = adv_val + D * lap * dt
        fn[0,:] = fn[1,:]
        fn[-1,:] = fn[-2,:]
        fn[:,0] = fn[:,1]
        fn[:,-1] = fn[:,-2]
        self.dopant_field, self.dopant_field_next = self.dopant_field_next, self.dopant_field

        if self.growth_progress > 0:
            interface_r = self.growth_progress
            band_inner = max(0.0, interface_r - INTERFACE_THICKNESS)
            inner_idx = int(band_inner*(n-1))
            inter_idx = int(interface_r*(n-1))
            if inter_idx > inner_idx:
                band = self.dopant_field[inner_idx:inter_idx, :]
                if band.size > 0:
                    incorporated = band.mean() * INCORPORATION_RATE * dt * band.size
                    self.total_incorporated_dopant += incorporated
                    band *= (1.0 - INCORPORATION_RATE * dt)
                    samples = min(10, inter_idx - inner_idx)
                    for _ in range(samples):
                        theta = random.uniform(0, 2*math.pi)
                        self.incorporated_positions.append(theta)
                    self.dopant_field[inner_idx:inter_idx, :] = band
        if len(self.incorporated_positions) > 5000:
            self.incorporated_positions = self.incorporated_positions[-3000:]

    def _dopant_uniformity_metric(self):
        if len(self.incorporated_positions) < 20:
            return 1.0
        angles = np.array(self.incorporated_positions[-1000:])
        ux = np.cos(angles).mean()
        uy = np.sin(angles).mean()
        r = math.sqrt(ux*ux + uy*uy)
        uniformity = 1.0 - r
        f = self.dopant_field
        mean = f.mean() + 1e-9
        std_ratio = f.std()/mean
        uniformity *= (1.0 - 0.4*std_ratio)
        return clamp(uniformity, 0.0, 1.0)

    def initializeGL(self):
        glClearColor(0.05, 0.05, 0.1, 1.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glLightfv(GL_LIGHT0, GL_POSITION, (5.0, 5.0, 10.0, 1.0))
        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, (0.8, 0.8, 0.8, 1.0))
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, w / h if h else 1, 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        glRotatef(self.view_angle, 0, 1, 0)
        self.draw_axes()
        self.draw_acoustic_levitators()
        if self.show_flow_lines:
            self.draw_flow_lines()
        self.draw_melt()
        if self.show_dopants:
            self.draw_melt_dopants()
        if self.growth_progress > 0:
            self.draw_crystal()

    def draw_axes(self):
        glDisable(GL_LIGHTING)
        glBegin(GL_LINES)
        glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(2, 0, 0)
        glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 2, 0)
        glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 2)
        glEnd()
        glEnable(GL_LIGHTING)

    def draw_acoustic_levitators(self):
        glDisable(GL_LIGHTING)
        t = self.melt_wobble * 0.5
        for idx, src in enumerate(self.acoustic_sources):
            glPushMatrix()
            glTranslatef(*src)
            glColor3f(0.65, 0.65, 0.7)
            glutSolidCube(0.22)
            glColor4f(0.5, 0.6, 1.0, 0.08 * self.acoustic_intensity)
            glBegin(GL_LINES)
            for ang in range(0, 360, 20):
                rad = math.radians(ang)
                x = math.cos(rad)
                z = math.sin(rad)
                amp = 0.6 + 0.2*math.sin(t + idx)
                glVertex3f(0, 0, 0)
                glVertex3f(x*amp, 0.0, z*amp)
            glEnd()
            glPopMatrix()
        glEnable(GL_LIGHTING)

    def draw_flow_lines(self):
        if not self.flow_lines or self.gravity < 0.5:
            return
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)
        time_factor = (math.sin(self.melt_wobble * 0.5) + 1) * 0.5
        for flow in self.flow_lines:
            x, y, z, direction, length, color = flow
            glBegin(GL_LINE_STRIP)
            animated_length = length * (0.5 + 0.5 * time_factor)
            segments = 12
            for i in range(segments + 1):
                t = i / segments
                t_anim = (t + time_factor) % 1.0
                curve_x = math.sin(t_anim * math.pi * 2) * 0.04
                curve_z = math.cos(t_anim * math.pi * 2) * 0.04
                px = x + direction[0] * animated_length * t_anim + curve_x
                py = y + direction[1] * animated_length * t_anim
                pz = z + direction[2] * animated_length * t_anim + curve_z
                alpha = 1.0 - t_anim
                glColor4f(color[0], color[1], color[2], color[3] * alpha)
                glVertex3f(px, py, pz)
            glEnd()
        glLineWidth(1.0)
        glEnable(GL_LIGHTING)

    def draw_melt(self):
        if self.growth_progress >= 1.0:
            return
        melt_size = 1.0 - self.growth_progress * 0.8
        wobble = math.sin(self.melt_wobble) * 0.1 if self.gravity < 0.5 else 0
        glPushMatrix()
        if self.gravity > 0.5:
            glScalef(1.05, 0.92, 1.05)
        temp_factor = (self.temperature - 500) / 1000
        if temp_factor > 0.7:
            glColor3f(1.0, 0.9, 0.3)
        elif temp_factor > 0.5:
            glColor3f(1.0, 0.6, 0.0)
        else:
            glColor3f(0.85, 0.25, 0.05)
        segments_theta = 40
        segments_phi = 28
        glDisable(GL_LIGHTING)
        glBegin(GL_TRIANGLES)
        t = self.melt_wobble
        for i in range(segments_phi):
            phi1 = math.pi * i / segments_phi
            phi2 = math.pi * (i+1) / segments_phi
            for j in range(segments_theta):
                theta1 = 2*math.pi * j / segments_theta
                theta2 = 2*math.pi * (j+1) / segments_theta
                tris = [
                    ((theta1, phi1), (theta2, phi1), (theta2, phi2)),
                    ((theta1, phi1), (theta2, phi2), (theta1, phi2))
                ]
                for (th_a, ph_a), (th_b, ph_b), (th_c, ph_c) in tris:
                    mod_a = self._surface_acoustic_modulation(th_a, ph_a, t)
                    mod_b = self._surface_acoustic_modulation(th_b, ph_b, t)
                    mod_c = self._surface_acoustic_modulation(th_c, ph_c, t)
                    ra = melt_size * self.melt_scale * (1.0 + wobble + mod_a)
                    rb = melt_size * self.melt_scale * (1.0 + wobble + mod_b)
                    rc = melt_size * self.melt_scale * (1.0 + wobble + mod_c)
                    ax = ra * math.sin(ph_a) * math.cos(th_a)
                    ay = ra * math.cos(ph_a)
                    az = ra * math.sin(ph_a) * math.sin(th_a)
                    bx = rb * math.sin(ph_b) * math.cos(th_b)
                    by = rb * math.cos(ph_b)
                    bz = rb * math.sin(ph_b) * math.sin(th_b)
                    cx = rc * math.sin(ph_c) * math.cos(th_c)
                    cy = rc * math.cos(ph_c)
                    cz = rc * math.sin(ph_c) * math.sin(th_c)
                    nx, ny, nz = compute_normal((ax, ay, az), (bx, by, bz), (cx, cy, cz))
                    glColor4f(1.0, 0.5 + 0.3*nx, 0.3 + 0.3*ny, 0.35)
                    glVertex3f(ax, ay, az); glVertex3f(bx, by, bz); glVertex3f(cx, cy, cz)
        glEnd()
        glEnable(GL_LIGHTING)
        glDisable(GL_LIGHTING)
        glColor4f(1.0, 0.8, 0.2, 0.15)
        glutSolidSphere(melt_size * 1.05 * (1.0 + wobble), 24, 24)
        glEnable(GL_LIGHTING)
        glPopMatrix()

    def draw_melt_dopants(self):
        if self.growth_progress >= 1.0:
            return
        glDisable(GL_LIGHTING)
        glColor4f(0.0, 0.9, 0.15, 0.75)
        n = DOPANT_N
        stride = max(1, n // 20)
        for i in range(0, n, stride):
            for j in range(0, n, stride):
                r = i/(n-1)
                z = j/(n-1)*DOPANT_HEIGHT - DOPANT_HEIGHT/2
                theta = 2*math.pi * ((i*37 + j*19) % 360)/360.0
                melt_radius = 1.0 - self.growth_progress*0.8
                if r <= melt_radius:
                    x = r * math.cos(theta)
                    y = z * 0.8
                    zc = r * math.sin(theta)
                    c = self.dopant_field[i,j]
                    size = 0.015 + 0.015 * (c)
                    glPushMatrix()
                    glTranslatef(x, y, zc)
                    glutSolidSphere(size, 6, 6)
                    glPopMatrix()
        glEnable(GL_LIGHTING)

    def draw_crystal(self):
        if self.growth_progress <= 0:
            return
        glPushMatrix()
        glTranslatef(0, -1.0 * self.growth_progress, 0)
        scale_factor = self.growth_progress
        glScalef(scale_factor, scale_factor, scale_factor)
        mat = get_material_data(self.material)
        crystal_color = mat["color"]
        brightness = 0.55 + 0.45 * self.crystal_quality
        glColor3f(crystal_color[0] * brightness,
                  crystal_color[1] * brightness,
                  crystal_color[2] * brightness)
        if self.gravity > 0.5:
            ir = 0.05 * (1.0 - self.crystal_quality)
            glScalef(1.0 + ir, 1.0 - ir * 2, 1.0 + ir)
        glBegin(GL_QUADS)
        for face in self.crystal_faces:
            if len(face) >= 3:
                v1 = self.crystal_vertices[face[0]]
                v2 = self.crystal_vertices[face[1]]
                v3 = self.crystal_vertices[face[2]]
                nx, ny, nz = compute_normal(v1, v2, v3)
                glNormal3f(nx, ny, nz)
                for idx in face:
                    glVertex3f(*self.crystal_vertices[idx])
        glEnd()
        if self.show_defects:
            self.draw_defects()
        if self.show_dopants:
            self.draw_crystal_dopants()
        glPopMatrix()

    def draw_defects(self):
        glDisable(GL_LIGHTING)
        for (x, y, z, size, defect_type) in self.defects:
            if defect_type == 'point':
                glColor4f(0.85, 0.05, 0.05, 0.7)
                glPushMatrix(); glTranslatef(x, y, z); glutSolidSphere(size, 8, 8); glPopMatrix()
            elif defect_type == 'line':
                glColor4f(1.0, 0.4, 0.05, 0.7)
                glPushMatrix(); glTranslatef(x, y, z)
                glBegin(GL_LINES); glVertex3f(0, -size * 3, 0); glVertex3f(0, size * 3, 0); glEnd()
                glPopMatrix()
            elif defect_type == 'cluster':
                glColor4f(0.75, 0.0, 0.75, 0.55)
                glPushMatrix(); glTranslatef(x, y, z)
                for _ in range(5):
                    dx = random.uniform(-1, 1) * size
                    dy = random.uniform(-1, 1) * size
                    dz = random.uniform(-1, 1) * size
                    glPushMatrix(); glTranslatef(dx, dy, dz); glutSolidSphere(size * 0.45, 6, 6); glPopMatrix()
                glPopMatrix()
        glEnable(GL_LIGHTING)

    def draw_crystal_dopants(self):
        glDisable(GL_LIGHTING)
        if self.incorporated_positions:
            glColor4f(0.1, 0.9, 0.3, 0.75)
            r = 0.45 * self.growth_progress
            for theta in self.incorporated_positions[-200:]:
                x = r * math.cos(theta)
                z = r * math.sin(theta)
                y = random.uniform(-0.4, 0.4) * self.growth_progress
                glPushMatrix(); glTranslatef(x, y, z); glutSolidSphere(0.02 * self.growth_progress, 6, 6); glPopMatrix()
        glEnable(GL_LIGHTING)

    def update_simulation(self):
        # PATCH: skip heavy work if hidden (prevents some driver errors when tabbing)
        if not self.isVisible():
            return
        self.view_angle += 0.1 * self.rotation_speed
        self.melt_wobble += 0.05
        if self.auto_grow:
            old = self.growth_progress
            self.growth_progress += self.auto_grow_speed * self.animation_direction
            self.growth_progress = clamp(self.growth_progress, 0.0, 1.0)
            if random.random() < 0.04:
                self.generate_defects()
                self.generate_dopants()
            if self.growth_progress != old:
                # Regenerate dopants to reflect new size snapshot
                self.generate_dopants()

        dt = min(MAX_DT, 0.016 + 0.002 * random.random())
        try:
            self._dopant_step(dt)
        except Exception as e:
            dprint("Dopant step error:", e)

        try:
            u_metric = self._dopant_uniformity_metric()
            self.displayed_uniformity = (1-UNIFORMITY_SMOOTH)*self.displayed_uniformity + UNIFORMITY_SMOOTH*u_metric
        except Exception as e:
            dprint("Uniformity metric error:", e)

        mat = get_material_data(self.material)
        ideal_temp = mat["melting_point"] - 50
        temp_diff = abs(self.temperature - ideal_temp)
        temp_score = clamp(1.0 - (temp_diff / 400.0), 0.0, 1.0)
        if self.gravity < 0.5:
            base_defect_factor = 1.0 - mat["defect_space"]
        else:
            base_defect_factor = 1.0 - mat["defect_earth"]
        defect_penalty = clamp(base_defect_factor, 0.0, 1.0)
        acoustic_score = 0.5 + 0.5 * self.acoustic_intensity
        effective_uniformity = clamp(self.displayed_uniformity + self.acoustic_intensity * ACOUSTIC_UNIFORMITY_BOOST, 0.0, 1.0)

        self.crystal_quality = (
            QUALITY_WEIGHT_DOPANT_UNIFORM * effective_uniformity +
            QUALITY_WEIGHT_TEMP * temp_score +
            QUALITY_WEIGHT_ACOUSTIC * acoustic_score +
            QUALITY_WEIGHT_DEFECT * defect_penalty
        )
        self.crystal_quality = clamp(self.crystal_quality, 0.0, 1.0)

        if self.gravity > 0.5 and random.random() < 0.01:
            self.generate_flow_lines()

        self.updateGL()


# Main Window

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle("NASA Quantum Material Crystal Growth Simulation (Enhanced Transport)")
        self.setMinimumSize(1400, 900)

        main_widget = QWidget()
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

        self.tab_widget = QTabWidget()
        self.simulation = CrystalGrowthSimulation()
        self.convection_view = ConvectionVisualizer(self.simulation)
        self.tab_widget.addTab(self.simulation, "3D Crystal Growth")
        self.tab_widget.addTab(self.convection_view, "Transport & Convection (2D)")
        main_layout.addWidget(self.tab_widget, 1)

        control_frame = QFrame()
        control_frame.setFrameShape(QFrame.StyledPanel)
        control_layout = QHBoxLayout()
        control_frame.setLayout(control_layout)
        main_layout.addWidget(control_frame, 0)

        # Environment
        env_group = QGroupBox("Environment")
        env_layout = QVBoxLayout()
        env_group.setLayout(env_layout)
        gravity_row = QHBoxLayout()
        self.earth_radio = QRadioButton("Earth (1g)")
        self.space_radio = QRadioButton("Space (µg)")
        self.space_radio.setChecked(True)
        gravity_row.addWidget(self.earth_radio)
        gravity_row.addWidget(self.space_radio)
        env_layout.addLayout(gravity_row)
        self.earth_radio.toggled.connect(self.update_gravity)

        temp_row = QHBoxLayout()
        temp_row.addWidget(QLabel("Temperature (°C):"))
        self.temp_slider = QSlider(Qt.Horizontal)
        self.temp_slider.setMinimum(500)
        self.temp_slider.setMaximum(1500)
        self.temp_slider.setValue(800)
        self.temp_label = QLabel("800°C")
        temp_row.addWidget(self.temp_slider)
        temp_row.addWidget(self.temp_label)
        self.temp_slider.valueChanged.connect(self.update_temperature)
        env_layout.addLayout(temp_row)

        acoustic_row = QHBoxLayout()
        acoustic_row.addWidget(QLabel("Acoustic Intensity:"))
        self.acoustic_slider = QSlider(Qt.Horizontal)
        self.acoustic_slider.setMinimum(0)
        self.acoustic_slider.setMaximum(100)
        self.acoustic_slider.setValue(50)
        self.acoustic_label = QLabel("50%")
        acoustic_row.addWidget(self.acoustic_slider)
        acoustic_row.addWidget(self.acoustic_label)
        self.acoustic_slider.valueChanged.connect(self.update_acoustic)
        env_layout.addLayout(acoustic_row)

        control_layout.addWidget(env_group)

        # Material
        material_group = QGroupBox("Material")
        material_layout = QVBoxLayout()
        material_group.setLayout(material_layout)
        self.material_combo = QComboBox()
        for m in QUANTUM_MATERIALS.keys():
            self.material_combo.addItem(m)
        self.material_combo.currentIndexChanged.connect(self.update_material)
        material_layout.addWidget(self.material_combo)
        self.material_info = QLabel("")
        self.material_info.setWordWrap(True)
        material_layout.addWidget(self.material_info)
        control_layout.addWidget(material_group)

        # Growth & Visualization
        growth_group = QGroupBox("Growth & Visualization")
        growth_layout = QVBoxLayout()
        growth_group.setLayout(growth_layout)
        growth_slider_row = QHBoxLayout()
        growth_slider_row.addWidget(QLabel("Growth Progress:"))
        self.growth_slider = QSlider(Qt.Horizontal)
        self.growth_slider.setMinimum(0)
        self.growth_slider.setMaximum(100)
        self.growth_slider.setValue(0)
        self.growth_value_label = QLabel("0%")
        growth_slider_row.addWidget(self.growth_slider)
        growth_slider_row.addWidget(self.growth_value_label)
        self.growth_slider.valueChanged.connect(self.update_growth)
        growth_layout.addLayout(growth_slider_row)
        self.auto_grow_btn = QPushButton("Start Auto Growth")
        self.auto_grow_btn.clicked.connect(self.toggle_auto_growth)
        growth_layout.addWidget(self.auto_grow_btn)
        rot_row = QHBoxLayout()
        rot_row.addWidget(QLabel("Rotation Speed:"))
        self.rotation_slider = QSlider(Qt.Horizontal)
        self.rotation_slider.setMinimum(0)
        self.rotation_slider.setMaximum(100)
        self.rotation_slider.setValue(50)
        self.rotation_label = QLabel("0.50")
        rot_row.addWidget(self.rotation_slider)
        rot_row.addWidget(self.rotation_label)
        self.rotation_slider.valueChanged.connect(self.update_rotation)
        growth_layout.addLayout(rot_row)
        vis_checks_layout = QVBoxLayout()
        self.flow_check = QCheckBox("Show Thermal Flow / Convection")
        self.flow_check.setChecked(True)
        self.flow_check.toggled.connect(self.update_visual_toggles)
        self.dopant_check = QCheckBox("Show Dopants")
        self.dopant_check.setChecked(True)
        self.dopant_check.toggled.connect(self.update_visual_toggles)
        self.defect_check = QCheckBox("Show Defects")
        self.defect_check.setChecked(True)
        self.defect_check.toggled.connect(self.update_visual_toggles)
        self.lattice_check = QCheckBox("Show Crystal Lattice Detail")
        self.lattice_check.setChecked(True)
        self.lattice_check.toggled.connect(self.toggle_lattice)
        for w in [self.flow_check, self.dopant_check, self.defect_check, self.lattice_check]:
            vis_checks_layout.addWidget(w)
        growth_layout.addLayout(vis_checks_layout)
        self.reset_btn = QPushButton("Reset")
        self.reset_btn.clicked.connect(self.reset_simulation)
        growth_layout.addWidget(self.reset_btn)
        control_layout.addWidget(growth_group)

        # Status
        status_group = QGroupBox("Status & Physics Explanation")
        status_layout = QVBoxLayout()
        status_group.setLayout(status_layout)
        self.quality_label = QLabel("Crystal Quality: 100.0%")
        self.defect_label = QLabel("Defects: 0")
        self.env_label = QLabel("Environment: Space (Microgravity)")
        self.dopant_label = QLabel("Dopant Uniformity: (calculating)")
        self.convection_label = QLabel("Convection: Diffusion Dominant")
        self.acoustic_label_status = QLabel("Acoustic Field: Intensity 50%")
        for w in [self.quality_label, self.defect_label, self.env_label,
                  self.dopant_label, self.convection_label, self.acoustic_label_status]:
            w.setWordWrap(True)
            status_layout.addWidget(w)
        self.explanation = QLabel(
            "<b>Earth (1g) Effects:</b><br>"
            "• Structured buoyant convection rolls advect dopants<br>"
            "• Dopant gradients & higher defect likelihood<br>"
            "• Acoustic field can partially suppress flow<br>"
            "<b>Space (Microgravity) Benefits:</b><br>"
            "• Transport governed mainly by diffusion<br>"
            "• More uniform dopant incorporation and fewer defects<br>"
            "• Acoustic levitation maintains containerless purity"
        )
        self.explanation.setWordWrap(True)
        status_layout.addWidget(self.explanation)
        control_layout.addWidget(status_group)

        self.status_timer = QTimer(self)
        self.status_timer.timeout.connect(self.refresh_status)
        self.status_timer.start(500)

        self.update_material(0)

    # UI Handlers
    def update_gravity(self):
        if self.earth_radio.isChecked():
            self.simulation.gravity = 1.0
            self.convection_view.simulation.gravity = 1.0
        else:
            self.simulation.gravity = 0.0
            self.convection_view.simulation.gravity = 0.0
        self.convection_view.reset_particles()
        self.simulation.generate_flow_lines()
        self.simulation.generate_defects()
        self.simulation.generate_dopants()

    def update_temperature(self):
        val = self.temp_slider.value()
        self.simulation.temperature = val
        self.temp_label.setText(f"{val}°C")

    def update_acoustic(self):
        val = self.acoustic_slider.value() / 100.0
        self.simulation.acoustic_intensity = val
        self.acoustic_label.setText(f"{int(val*100)}%")

    def update_material(self, idx):
        material = self.material_combo.currentText()
        self.simulation.material = material
        self.simulation.generate_crystal_structure()
        self.simulation.generate_dopants()
        data = get_material_data(material)
        info = (f"Melting Point: {data['melting_point']}°C\n"
                f"Structure: {data['crystal_structure']}\n"
                f"Nominal Growth Rate: {data['growth_rate']} mm/h\n"
                f"Earth Defect Prob: {data['defect_earth']*100:.1f}%\n"
                f"Space Defect Prob: {data['defect_space']*100:.1f}%")
        self.material_info.setText(info)

    def update_growth(self):
        val = self.growth_slider.value() / 100.0
        self.simulation.growth_progress = clamp(val,0.0,1.0)
        self.growth_value_label.setText(f"{int(self.simulation.growth_progress*100)}%")
        self.simulation.generate_defects()
        self.simulation.generate_dopants()

    def update_rotation(self):
        val = self.rotation_slider.value() / 100.0
        self.simulation.rotation_speed = val
        self.rotation_label.setText(f"{val:.2f}")

    def update_visual_toggles(self):
        self.simulation.show_flow_lines = self.flow_check.isChecked()
        self.simulation.show_dopants = self.dopant_check.isChecked()
        self.simulation.show_defects = self.defect_check.isChecked()
        self.convection_view.show_convection = self.flow_check.isChecked()
        self.convection_view.show_dopants = self.dopant_check.isChecked()

    def toggle_lattice(self):
        self.simulation.show_crystal_lattice = self.lattice_check.isChecked()
        self.simulation.generate_crystal_structure()

    def toggle_auto_growth(self):
        self.simulation.auto_grow = not self.simulation.auto_grow
        self.auto_grow_btn.setText("Pause Auto Growth" if self.simulation.auto_grow else "Start Auto Growth")

    def reset_simulation(self):
        self.simulation.growth_progress = 0.0
        self.growth_slider.setValue(0)
        self.simulation.generate_defects()
        self.simulation.generate_dopants()
        self.simulation.auto_grow = False
        self.auto_grow_btn.setText("Start Auto Growth")
        self.simulation.crystal_quality = 1.0
        self.simulation._init_dopant_field()
        self.simulation.incorporated_positions.clear()
        self.simulation.total_incorporated_dopant = 0.0

    def refresh_status(self):
        quality = self.simulation.crystal_quality * 100
        self.quality_label.setText(f"Crystal Quality: {quality:.1f}%")
        defect_count = len(self.simulation.defects)
        self.defect_label.setText(f"Defects: {defect_count}")
        env_text = "Space (Microgravity)" if self.simulation.gravity < 0.5 else "Earth (1g)"
        self.env_label.setText(f"Environment: {env_text}")
        uniformity = self.simulation.displayed_uniformity * 100
        self.dopant_label.setText(f"Dopant Uniformity: {uniformity:.1f}%")
        if self.simulation.gravity < 0.5:
            self.convection_label.setText("Convection: Diffusion Dominant")
        else:
            self.convection_label.setText("Convection: Structured Buoyant Rolls")
        self.acoustic_label_status.setText(f"Acoustic Field: Intensity {int(self.simulation.acoustic_intensity*100)}%")

# Entry
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

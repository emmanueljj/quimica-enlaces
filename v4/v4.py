import os
import re
import math
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import numpy as np
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ==========================================
# BACKEND 1: Generador de Coordenadas
# ==========================================
class SeparadorMoleculas:
    def __init__(self, archivo_gjf):
        self.archivo_gjf = archivo_gjf
        self.encabezado, self.titulo, self.pie = [], [], []
        self.carga_mult = ""
        self.coordenadas = []
        self._parsear_archivo()

    def _parsear_archivo(self):
        if not os.path.exists(self.archivo_gjf): return
        with open(self.archivo_gjf, 'r', encoding='utf-8') as f:
            lineas = f.readlines()

        seccion = 0
        for linea in lineas:
            linea_limpia = linea.strip()
            if linea_limpia == "" and seccion < 3:
                seccion += 1
                continue
            if seccion == 0: self.encabezado.append(linea)
            elif seccion == 1: self.titulo.append(linea)
            elif seccion == 2:
                if not self.carga_mult: self.carga_mult = linea
                else:
                    partes = linea_limpia.split()
                    if len(partes) >= 4 and not partes[0].isdigit():
                        self.coordenadas.append({
                            'elemento': partes[0], 'x': float(partes[1]), 
                            'y': float(partes[2]), 'z': float(partes[3])
                        })
                    else:
                        seccion = 3
                        self.pie.append(linea)
            elif seccion >= 3: self.pie.append(linea)

    def get_distancia_actual(self, atomos_mol1):
        if atomos_mol1 >= len(self.coordenadas) or atomos_mol1 <= 0: return 0.0
        mol1 = self.coordenadas[:atomos_mol1]
        c1 = [sum(a[e] for a in mol1)/len(mol1) for e in ['x','y','z']]
        mol2 = self.coordenadas[atomos_mol1:]
        c2 = [sum(a[e] for a in mol2)/len(mol2) for e in ['x','y','z']]
        return math.sqrt(sum((c2[i]-c1[i])**2 for i in range(3)))

    def generar_separacion(self, atomos_mol1, dist_ini, dist_fin, avance, directorio_salida, nombre_base, callback=None):
        if not os.path.exists(directorio_salida): os.makedirs(directorio_salida)
        mol1, mol2 = self.coordenadas[:atomos_mol1], self.coordenadas[atomos_mol1:]
        c1 = [sum(a[e] for a in mol1)/len(mol1) for e in ['x','y','z']]
        c2 = [sum(a[e] for a in mol2)/len(mol2) for e in ['x','y','z']]
        
        v = [c2[i] - c1[i] for i in range(3)]
        mag_orig = math.sqrt(sum(coord**2 for coord in v))
        if mag_orig == 0: raise ValueError("Centroides coinciden. No se puede definir dirección.")
        u = [coord / mag_orig for coord in v]

        targets = []
        d = dist_ini
        avance = abs(avance)
        if dist_ini <= dist_fin:
            while d <= dist_fin + 1e-5:
                targets.append(d); d += avance
        else:
            while d >= dist_fin - 1e-5:
                targets.append(d); d -= avance

        total = len(targets)
        for paso, dist_obj in enumerate(targets):
            ruta = os.path.join(directorio_salida, f"{nombre_base}_{paso}.gjf")
            with open(ruta, 'w', encoding='utf-8') as f:
                f.writelines(self.encabezado); f.write("\n")
                tit = self.titulo[0].strip() if self.titulo else "Generado"
                f.write(f"{tit} - Distancia {dist_obj:.3f}A - Paso {paso}\n\n")
                f.write(self.carga_mult)
                
                despl = dist_obj - mag_orig
                for i, at in enumerate(self.coordenadas):
                    if i < atomos_mol1:
                        f.write(f" {at['elemento']:<2} {at['x']:>14.8f} {at['y']:>14.8f} {at['z']:>14.8f}\n")
                    else:
                        x, y, z = at['x'] + u[0]*despl, at['y'] + u[1]*despl, at['z'] + u[2]*despl
                        f.write(f" {at['elemento']:<2} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")
                f.write("\n")
                if self.pie: f.writelines(self.pie); f.write("\n")
            if callback: callback(paso + 1, total)
        return total

# ==========================================
# BACKEND 2: Analizador de Energías
# ==========================================
class ExtractorGaussian:
    @staticmethod
    def extraer_energia(ruta_archivo):
        patron = re.compile(r"SCF Done:\s+E\([A-Za-z0-9]+\)\s*=\s*(-?\d+\.\d+)")
        try:
            with open(ruta_archivo, 'r', encoding='utf-8', errors='ignore') as f:
                for linea in f:
                    match = patron.search(linea)
                    if match: return float(match.group(1))
        except: pass
        return None

    @staticmethod
    def extraer_distancia_centroides(ruta_archivo, atomos_mol1):
        try:
            with open(ruta_archivo, 'r', encoding='utf-8', errors='ignore') as f: lineas = f.readlines()
        except: return None
        
        idx_inicio = -1
        for i in range(len(lineas)-1, -1, -1):
            if "Standard orientation:" in lineas[i] or "Input orientation:" in lineas[i]:
                idx_inicio = i; break
        if idx_inicio == -1: return None
            
        idx, coordenadas = idx_inicio + 5, []
        while idx < len(lineas) and "---------------------------------------------------------------------" not in lineas[idx]:
            partes = lineas[idx].split()
            if len(partes) >= 6: coordenadas.append([float(partes[3]), float(partes[4]), float(partes[5])])
            idx += 1
            
        if len(coordenadas) <= atomos_mol1: return None
        mol1, mol2 = coordenadas[:atomos_mol1], coordenadas[atomos_mol1:]
        c1 = [sum(a[i] for a in mol1)/len(mol1) for i in range(3)]
        c2 = [sum(a[i] for a in mol2)/len(mol2) for i in range(3)]
        return math.sqrt(sum((c2[i] - c1[i])**2 for i in range(3)))

    @staticmethod
    def extraer_scan_interno(ruta_archivo, atomos_mol1):
        """Lee un archivo .out gigante y extrae todos los puntos optimizados de un escaneo."""
        try:
            with open(ruta_archivo, 'r', encoding='utf-8', errors='ignore') as f:
                lineas = f.readlines()
        except Exception: return []
            
        puntos = []
        is_relaxed = any("Optimization completed." in l for l in lineas)
        
        dist_actual = None
        e_actual = None
        
        idx = 0
        while idx < len(lineas):
            linea = lineas[idx]
            
            # 1. Buscar Coordenadas
            if "Standard orientation:" in linea or "Input orientation:" in linea:
                idx += 5
                coordenadas = []
                while idx < len(lineas) and "---------------------------------------------------------------------" not in lineas[idx]:
                    partes = lineas[idx].split()
                    if len(partes) >= 6:
                        coordenadas.append([float(partes[3]), float(partes[4]), float(partes[5])])
                    idx += 1
                
                if len(coordenadas) >= atomos_mol1:
                    mol1 = coordenadas[:atomos_mol1]
                    mol2 = coordenadas[atomos_mol1:]
                    c1 = [sum(a[i] for a in mol1)/len(mol1) for i in range(3)]
                    c2 = [sum(a[i] for a in mol2)/len(mol2) for i in range(3)]
                    dist_actual = math.sqrt(sum((c2[i] - c1[i])**2 for i in range(3)))
                    
            # 2. Buscar Energía
            elif "SCF Done:" in linea:
                match = re.search(r"SCF Done:\s+E\([A-Za-z0-9]+\)\s*=\s*(-?\d+\.\d+)", linea)
                if match:
                    e_actual = float(match.group(1))
                    if not is_relaxed and dist_actual is not None:
                        puntos.append((dist_actual, e_actual))
                        dist_actual = None 
                        
            # 3. Confirmar Optimización (Solo en Scan Relaxed)
            elif is_relaxed and "Optimization completed." in linea:
                if dist_actual is not None and e_actual is not None:
                    puntos.append((dist_actual, e_actual))
                    dist_actual = None
                    e_actual = None
                    
            idx += 1
            
        return puntos

# ==========================================
# FRONTEND: Interfaz de Usuario Limpia
# ==========================================
class GaussianSuiteApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Suite de Análisis Molecular (Gaussian)")
        self.root.geometry("1050x700")
        
        self._configurar_estilos()
        
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.tab_gen = ttk.Frame(self.notebook)
        self.tab_graf = ttk.Frame(self.notebook)
        
        self.notebook.add(self.tab_gen, text="  🛠️ 1. Generador de Coordenadas (.gjf)  ")
        self.notebook.add(self.tab_graf, text="  📈 2. Análisis y Gráficas (.out)  ")
        
        self._init_tab_generador()
        self._init_tab_graficador()

    def _configurar_estilos(self):
        self.style = ttk.Style()
        self.style.theme_use('clam')
        # Configuraciones mínimas para resaltar texto o botones
        self.style.configure('Action.TButton', font=('Helvetica', 11, 'bold'))
        self.style.configure('TLabelframe.Label', font=('Helvetica', 11, 'bold'))

    # ---------------------------------------------------------
    # UI PESTAÑA 1: GENERADOR
    # ---------------------------------------------------------
    def _init_tab_generador(self):
        self.gen_ruta = tk.StringVar()
        self.gen_separador = None
        
        main_frame = ttk.Frame(self.tab_gen, padding=20)
        main_frame.pack(fill="both", expand=True)
        
        ttk.Label(main_frame, text="Configuración del Escaneo de Superficie", font=('Helvetica', 16, 'bold')).pack(anchor='w', pady=(0, 15))
        
        frm_archivo = ttk.LabelFrame(main_frame, text=" 📂 Archivo Plantilla (.gjf) ")
        frm_archivo.pack(fill="x", pady=5, ipady=10)
        
        frm_arch_in = ttk.Frame(frm_archivo)
        frm_arch_in.pack(fill='x', padx=15, pady=5)
        ttk.Entry(frm_arch_in, textvariable=self.gen_ruta, state="readonly").pack(side="left", fill="x", expand=True, padx=(0, 10))
        ttk.Button(frm_arch_in, text="Buscar Archivo...", command=self._gen_seleccionar).pack(side="right")
        
        self.gen_lbl_info = ttk.Label(frm_archivo, text="No se ha cargado ningún archivo.", font=('Helvetica', 9, 'italic'))
        self.gen_lbl_info.pack(anchor='w', padx=15)

        frm_param = ttk.LabelFrame(main_frame, text=" ⚙️ Parámetros Físicos ")
        frm_param.pack(fill="x", pady=15, ipady=10)
        
        grid_frame = ttk.Frame(frm_param)
        grid_frame.pack(padx=15, pady=5)
        
        pad_opts = {'padx': 10, 'pady': 8}
        
        ttk.Label(grid_frame, text="Átomos en Molécula 1:").grid(row=0, column=0, sticky="e", **pad_opts)
        self.gen_spin_m1 = ttk.Spinbox(grid_frame, from_=1, to=999, width=10, command=self._gen_act_dist)
        self.gen_spin_m1.bind('<KeyRelease>', lambda e: self._gen_act_dist())
        self.gen_spin_m1.grid(row=0, column=1, sticky="w", **pad_opts)
        
        ttk.Label(grid_frame, text="Distancia Inicial (Å):").grid(row=1, column=0, sticky="e", **pad_opts)
        self.gen_ent_ini = ttk.Entry(grid_frame, width=12); self.gen_ent_ini.insert(0, "1.5")
        self.gen_ent_ini.grid(row=1, column=1, sticky="w", **pad_opts)

        ttk.Label(grid_frame, text="Distancia Final (Å):").grid(row=1, column=2, sticky="e", **pad_opts)
        self.gen_ent_fin = ttk.Entry(grid_frame, width=12); self.gen_ent_fin.insert(0, "6.0")
        self.gen_ent_fin.grid(row=1, column=3, sticky="w", **pad_opts)

        ttk.Label(grid_frame, text="Avance (Å):").grid(row=1, column=4, sticky="e", **pad_opts)
        self.gen_ent_av = ttk.Entry(grid_frame, width=12); self.gen_ent_av.insert(0, "0.25")
        self.gen_ent_av.grid(row=1, column=5, sticky="w", **pad_opts)

        frm_out = ttk.LabelFrame(main_frame, text=" 📁 Opciones de Salida ")
        frm_out.pack(fill="x", pady=5, ipady=10)
        ttk.Label(frm_out, text="Prefijo archivos:").pack(side="left", padx=(15, 5))
        self.gen_ent_nom = ttk.Entry(frm_out, width=30); self.gen_ent_nom.insert(0, "scan_pes")
        self.gen_ent_nom.pack(side="left", padx=5)

        self.gen_prog = ttk.Progressbar(main_frame, orient="horizontal", mode="determinate")
        self.gen_prog.pack(fill="x", pady=(25, 10))
        
        self.gen_btn_run = ttk.Button(main_frame, text="GENERAR SISTEMAS GJF", style='Action.TButton', command=self._gen_ejecutar)
        self.gen_btn_run.pack(pady=5, ipadx=20)

    def _gen_seleccionar(self):
        ruta = filedialog.askopenfilename(title="Seleccionar .gjf", filetypes=[("Gaussian", "*.gjf")])
        if ruta:
            self.gen_ruta.set(ruta)
            try:
                self.gen_separador = SeparadorMoleculas(ruta)
                n = len(self.gen_separador.coordenadas)
                self.gen_spin_m1.config(to=n-1)
                if self.gen_spin_m1.get() == "1" and n > 2: self.gen_spin_m1.set(int(n/2))
                self._gen_act_dist()
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo leer:\n{e}")

    def _gen_act_dist(self):
        if self.gen_separador:
            try:
                m1 = int(self.gen_spin_m1.get())
                d = self.gen_separador.get_distancia_actual(m1)
                txt = f"✔️ {len(self.gen_separador.coordenadas)} Átomos | Distancia en archivo: {d:.4f} Å" if d > 0 else "Distancia inválida."
                self.gen_lbl_info.config(text=txt, foreground='green' if d>0 else 'red')
            except ValueError: pass

    def _gen_ejecutar(self):
        if not self.gen_separador: return messagebox.showwarning("Datos", "Selecciona un archivo válido.")
        try:
            m1 = int(self.gen_spin_m1.get())
            d_i, d_f, av = float(self.gen_ent_ini.get()), float(self.gen_ent_fin.get()), float(self.gen_ent_av.get())
            nom = self.gen_ent_nom.get().strip()
            if av <= 0 or d_i == d_f: raise ValueError("Avance debe ser > 0 y distancias distintas.")
            
            out_dir = os.path.join(os.path.dirname(self.gen_ruta.get()), f"scan_{nom}")
            self.gen_btn_run.config(state="disabled", text="Procesando...")
            
            total = self.gen_separador.generar_separacion(m1, d_i, d_f, av, out_dir, nom, 
                                                          lambda a, t: self._update_prog(self.gen_prog, a, t))
            if total > 0: messagebox.showinfo("Éxito", f"{total} archivos creados en:\n{out_dir}")
        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            self.gen_btn_run.config(state="normal", text="GENERAR SISTEMAS GJF")
            self.gen_prog['value'] = 0

    def _update_prog(self, pbar, act, tot):
        pbar['value'] = (act/tot)*100
        self.root.update_idletasks()

    # ---------------------------------------------------------
    # UI PESTAÑA 2: GRAFICADOR
    # ---------------------------------------------------------
    def _init_tab_graficador(self):
        self.graf_arch_sis, self.graf_arch_m1, self.graf_arch_m2 = [], [], []
        self.graf_datos = []
        self.graf_modo = 'multiple'
        
        main_paned = ttk.PanedWindow(self.tab_graf, orient=tk.HORIZONTAL)
        main_paned.pack(fill="both", expand=True, padx=10, pady=10)
        
        # --- PANEL IZQUIERDO ---
        panel_izq = ttk.Frame(main_paned, width=330)
        main_paned.add(panel_izq, weight=0)
        
        ttk.Label(panel_izq, text="Carga de Resultados (.out)", font=('Helvetica', 13, 'bold')).pack(anchor='w', pady=(0, 10))
        
        frm_sis = ttk.LabelFrame(panel_izq, text=" 1. Sistemas (Complejo) ")
        frm_sis.pack(fill="x", pady=5)
        
        btn_frm = ttk.Frame(frm_sis)
        btn_frm.pack(fill='x', padx=10, pady=5)
        ttk.Button(btn_frm, text="📄 1 Archivo (Scan)", command=self._graf_load_sis_single).pack(side='left', expand=True, fill='x', padx=(0,2))
        ttk.Button(btn_frm, text="📂 Puntos Sueltos", command=self._graf_load_sis_multiple).pack(side='left', expand=True, fill='x', padx=(2,0))
        
        self.lbl_sis = ttk.Label(frm_sis, text="0 archivos", font=('Helvetica', 9))
        self.lbl_sis.pack(pady=(0,5))

        def crear_bloque_carga(padre, texto, cmd):
            frm = ttk.LabelFrame(padre, text=texto)
            frm.pack(fill="x", pady=5)
            btn = ttk.Button(frm, text="Cargar Archivo(s)", command=cmd)
            btn.pack(pady=5, fill='x', padx=10)
            lbl = ttk.Label(frm, text="0 archivos", font=('Helvetica', 9))
            lbl.pack(pady=(0,5))
            return lbl

        self.lbl_m1 = crear_bloque_carga(panel_izq, " 2. Ref. Molécula 1 Aislada ", self._graf_load_m1)
        self.lbl_m2 = crear_bloque_carga(panel_izq, " 3. Ref. Molécula 2 Aislada ", self._graf_load_m2)

        frm_param = ttk.LabelFrame(panel_izq, text=" Parámetros Analíticos ")
        frm_param.pack(fill="x", pady=10)
        
        frm_spin = ttk.Frame(frm_param)
        frm_spin.pack(pady=5)
        ttk.Label(frm_spin, text="Átomos Mol 1:").pack(side="left", padx=5)
        self.graf_spin_m1 = ttk.Spinbox(frm_spin, from_=1, to=999, width=5); self.graf_spin_m1.set(2)
        self.graf_spin_m1.pack(side="left")
        
        self.graf_var_lj = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm_param, text="Ajuste Lennard-Jones (L-J)", variable=self.graf_var_lj).pack(pady=10)

        self.graf_btn_calc = ttk.Button(panel_izq, text="CALCULAR Y GRAFICAR", style='Action.TButton', command=self._graf_procesar)
        self.graf_btn_calc.pack(fill="x", pady=15, ipady=4)

        self.graf_btn_save = ttk.Button(panel_izq, text="💾 Guardar Imagen (.png)", command=self._graf_guardar, state='disabled')
        self.graf_btn_save.pack(fill="x")

        # --- PANEL DERECHO (Gráfica) ---
        panel_der = ttk.Frame(main_paned)
        main_paned.add(panel_der, weight=1)
        
        self.fig, self.ax = plt.subplots(figsize=(6, 5), dpi=100)
        self.ax.set_title("Esperando datos...", fontsize=12, color='#666')
        self.ax.grid(True, linestyle='--', alpha=0.5)
        self.ax.set_axis_off() 
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=panel_der)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=5)

    # Lógica Carga Archivos
    def _graf_load_sis_single(self):
        arch = filedialog.askopenfilename(title="Seleccionar Escaneo Único (.out)", filetypes=[("Archivos OUT", "*.out *.log")])
        if arch:
            self.graf_arch_sis = [arch]
            self.graf_modo = 'single'
            self.lbl_sis.config(text=f"1 Archivo Escaneo Único", foreground='blue')
            
    def _graf_load_sis_multiple(self):
        archs = filedialog.askopenfilenames(title="Seleccionar Puntos Sueltos", filetypes=[("Archivos OUT", "*.out *.log")])
        if archs:
            self.graf_arch_sis = list(archs)
            self.graf_modo = 'multiple'
            self.lbl_sis.config(text=f"{len(archs)} archivos sueltos", foreground='blue')
            
    def _graf_load_m1(self):
        archs = filedialog.askopenfilenames(filetypes=[("Archivos OUT", "*.out *.log")])
        if archs: self.graf_arch_m1 = list(archs); self.lbl_m1.config(text=f"{len(archs)} archivos", foreground='blue')
    def _graf_load_m2(self):
        archs = filedialog.askopenfilenames(filetypes=[("Archivos OUT", "*.out *.log")])
        if archs: self.graf_arch_m2 = list(archs); self.lbl_m2.config(text=f"{len(archs)} archivos", foreground='blue')

    def _graf_guardar(self):
        if not self.graf_datos: return
        ruta = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png")])
        if ruta:
            try:
                self.fig.savefig(ruta, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Guardado", "Imagen guardada exitosamente.")
            except Exception as e: messagebox.showerror("Error", str(e))

    def _graf_procesar(self):
        if not self.graf_arch_sis or not self.graf_arch_m1 or not self.graf_arch_m2:
            return messagebox.showwarning("Datos", "Carga los archivos de los 3 bloques.")
        try: m1_at = int(self.graf_spin_m1.get())
        except: return messagebox.showerror("Error", "Átomos inválidos.")

        self.graf_btn_calc.config(state='disabled', text='Procesando...')
        self.root.update()

        try:
            e_m1 = [ExtractorGaussian.extraer_energia(f) for f in self.graf_arch_m1]
            e_m2 = [ExtractorGaussian.extraer_energia(f) for f in self.graf_arch_m2]
            if None in e_m1 or None in e_m2: raise ValueError("Falta energía (SCF Done) en archivos de referencia aislada.")
            
            e1_uni = e_m1[0] if len(e_m1)==1 else None
            e2_uni = e_m2[0] if len(e_m2)==1 else None
            
            self.graf_datos = []
            HARTREE_TO_KCAL = 627.509
            
            if self.graf_modo == 'single':
                if e1_uni is None or e2_uni is None:
                    raise ValueError("Para graficar un Escaneo Único, necesitas proveer exactamente 1 archivo de referencia para M1 y 1 para M2.")
                    
                puntos_escaneo = ExtractorGaussian.extraer_scan_interno(self.graf_arch_sis[0], m1_at)
                if not puntos_escaneo:
                    raise ValueError("No se encontraron geometrías o energías válidas en el archivo del escaneo.")
                
                for dist, e_s in puntos_escaneo:
                    self.graf_datos.append((dist, (e_s - (e1_uni + e2_uni)) * HARTREE_TO_KCAL))
            else:
                for i, r_sis in enumerate(self.graf_arch_sis):
                    e_s = ExtractorGaussian.extraer_energia(r_sis)
                    d = ExtractorGaussian.extraer_distancia_centroides(r_sis, m1_at)
                    if e_s is None or d is None: continue
                    e1 = e1_uni if e1_uni is not None else e_m1[i]
                    e2 = e2_uni if e2_uni is not None else e_m2[i]
                    self.graf_datos.append((d, (e_s - (e1 + e2)) * HARTREE_TO_KCAL)) 

            if not self.graf_datos: raise ValueError("No se extrajeron datos válidos.")
            
            self.graf_datos.sort(key=lambda x: x[0])
            dists, ener = [d[0] for d in self.graf_datos], [d[1] for d in self.graf_datos]

            self.ax.clear()
            self.ax.set_axis_on()
            
            self.ax.plot(dists, ener, marker='o', ls='-', c='#D32F2F', lw=2, label="DFT (Gaussian)")

            if self.graf_var_lj.get():
                try:
                    def lj_despl(r, A, B, C): return (A/r**12) - (B/r**6) + C
                    popt, _ = curve_fit(lj_despl, np.array(dists), np.array(ener), maxfev=100000, bounds=([0,0,-np.inf], [np.inf, np.inf, np.inf]))
                    r_sm = np.linspace(min(dists)*0.95, max(dists)*1.05, 200)
                    self.ax.plot(r_sm, lj_despl(r_sm, *popt), ls='--', c='#1976D2', lw=2.5, label="Ajuste L-J")
                    print(f"Parámetros LJ: A={popt[0]:.2e}, B={popt[1]:.2e}, C={popt[2]:.2f}")
                except Exception as e: print(f"Fallo ajuste LJ: {e}")

            self.ax.set_title("Curva de Energía de Interacción Molecular", fontsize=14, fontweight='bold', pad=15)
            self.ax.set_xlabel("Distancia entre Centroides (Å)", fontsize=11)
            self.ax.set_ylabel("ΔE [kcal/mol]", fontsize=11)
            self.ax.grid(True, ls='--', alpha=0.6)
            self.ax.axhline(0, color='black', ls=':', alpha=0.7)
            self.ax.legend()
            
            y_min, y_max = min(ener), max(ener)
            margen = (y_max - y_min)*0.05 if y_max != y_min else 1.0
            self.ax.set_ylim(y_min - margen, y_max + margen)
            self.canvas.draw()
            
            self.graf_btn_save.config(state='normal')
            messagebox.showinfo("Éxito", f"Se graficaron {len(self.graf_datos)} puntos exitosamente.")
            
        except Exception as e: messagebox.showerror("Error de Análisis", str(e))
        finally: self.graf_btn_calc.config(state='normal', text='CALCULAR Y GRAFICAR')

if __name__ == "__main__":
    root = tk.Tk()
    app = GaussianSuiteApp(root)
    root.mainloop()
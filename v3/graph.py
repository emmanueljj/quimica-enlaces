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

class ExtractorGaussian:
    @staticmethod
    def extraer_energia(ruta_archivo):
        energia = None
        patron_energia = re.compile(r"SCF Done:\s+E\([A-Za-z0-9]+\)\s*=\s*(-?\d+\.\d+)")
        try:
            with open(ruta_archivo, 'r', encoding='utf-8', errors='ignore') as f:
                for linea in f:
                    match = patron_energia.search(linea)
                    if match:
                        energia = float(match.group(1))
        except Exception as e:
            print(f"Error leyendo {ruta_archivo}: {e}")
        return energia

    @staticmethod
    def extraer_distancia_centroides(ruta_archivo, atomos_mol1):
        lineas = []
        try:
            with open(ruta_archivo, 'r', encoding='utf-8', errors='ignore') as f:
                lineas = f.readlines()
        except Exception:
            return None
            
        idx_inicio = -1
        for i in range(len(lineas)-1, -1, -1):
            if "Standard orientation:" in lineas[i] or "Input orientation:" in lineas[i]:
                idx_inicio = i
                break
                
        if idx_inicio == -1: return None
            
        idx = idx_inicio + 5
        coordenadas = []
        while idx < len(lineas) and "---------------------------------------------------------------------" not in lineas[idx]:
            partes = lineas[idx].split()
            if len(partes) >= 6:
                coordenadas.append([float(partes[3]), float(partes[4]), float(partes[5])])
            idx += 1
            
        if len(coordenadas) <= atomos_mol1: return None
            
        mol1 = coordenadas[:atomos_mol1]
        mol2 = coordenadas[atomos_mol1:]
        
        c1 = [sum(a[i] for a in mol1) / len(mol1) for i in range(3)]
        c2 = [sum(a[i] for a in mol2) / len(mol2) for i in range(3)]
        
        return math.sqrt(sum((c2[i] - c1[i])**2 for i in range(3)))

class AppGraficador:
    def __init__(self, root):
        self.root = root
        self.root.title("Análisis de Interacción y Lennard-Jones")
        self.root.geometry("950x700")
        
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        self.archivos_sis = []
        self.archivos_m1 = []
        self.archivos_m2 = []
        self.datos_grafica = [] 
        
        self.var_sis = tk.StringVar(value="0 archivos seleccionados")
        self.var_m1 = tk.StringVar(value="0 archivos seleccionados")
        self.var_m2 = tk.StringVar(value="0 archivos seleccionados")
        self.var_mostrar_lj = tk.BooleanVar(value=True) 
        
        self.crear_interfaz()

    def crear_interfaz(self):
        panel_ctrl = ttk.Frame(self.root, width=350)
        panel_ctrl.pack(side="left", fill="y", padx=15, pady=15)
        
        ttk.Label(panel_ctrl, text="Extracción de Datos", font=("Helvetica", 14, "bold")).pack(pady=(0, 10))
        
        frame_sis = ttk.LabelFrame(panel_ctrl, text=" 1. Sistemas Completos ")
        frame_sis.pack(fill="x", pady=5)
        ttk.Button(frame_sis, text="Cargar Puntos de Escaneo", command=self.cargar_sis).pack(pady=5)
        ttk.Label(frame_sis, textvariable=self.var_sis, foreground="blue").pack(pady=2)
        
        frame_m1 = ttk.LabelFrame(panel_ctrl, text=" 2. Molécula 1 Aislada ")
        frame_m1.pack(fill="x", pady=5)
        ttk.Button(frame_m1, text="Cargar Referencia M1", command=self.cargar_m1).pack(pady=5)
        ttk.Label(frame_m1, textvariable=self.var_m1, foreground="blue").pack(pady=2)
        
        frame_m2 = ttk.LabelFrame(panel_ctrl, text=" 3. Molécula 2 Aislada ")
        frame_m2.pack(fill="x", pady=5)
        ttk.Button(frame_m2, text="Cargar Referencia M2", command=self.cargar_m2).pack(pady=5)
        ttk.Label(frame_m2, textvariable=self.var_m2, foreground="blue").pack(pady=2)
        
        frame_param = ttk.LabelFrame(panel_ctrl, text=" Parámetros y Opciones ")
        frame_param.pack(fill="x", pady=10)
        
        frame_spin = ttk.Frame(frame_param)
        frame_spin.pack(pady=5)
        ttk.Label(frame_spin, text="Átomos Mol 1:").pack(side="left", padx=5)
        self.spin_m1 = ttk.Spinbox(frame_spin, from_=1, to=999, width=5)
        self.spin_m1.set(2)
        self.spin_m1.pack(side="left")
        
        ttk.Checkbutton(frame_param, text="Ajustar Potencial de Lennard-Jones", variable=self.var_mostrar_lj).pack(pady=5, anchor="w", padx=10)
        
        self.btn_calc = ttk.Button(panel_ctrl, text="CALCULAR Y GRAFICAR", command=self.procesar_y_graficar)
        self.style.configure('Action.TButton', font=('Helvetica', 11, 'bold'), background='#4CAF50', foreground='white')
        self.btn_calc.configure(style='Action.TButton')
        self.btn_calc.pack(fill="x", pady=15, ipady=5)

        self.btn_guardar = ttk.Button(panel_ctrl, text="💾 GUARDAR GRÁFICA (.PNG)", command=self.guardar_grafica, state="disabled")
        self.btn_guardar.pack(fill="x", pady=5, ipady=3)

        # --- PANEL DERECHO (Gráfica) ---
        self.panel_grafica = ttk.Frame(self.root)
        self.panel_grafica.pack(side="right", fill="both", expand=True, padx=10, pady=10)
        
        self.fig, self.ax = plt.subplots(figsize=(6, 5), dpi=100)
        self.ax.set_title("Curva de Energía de Interacción", fontsize=12)
        self.ax.set_xlabel("Distancia (Å)")
        self.ax.set_ylabel("Energía de Interacción (kcal/mol)")
        self.ax.grid(True, linestyle='--', alpha=0.7)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.panel_grafica)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def cargar_sis(self):
        archivos = filedialog.askopenfilenames(title="Sistemas Completos", filetypes=[("Archivos OUT", "*.out *.log")])
        if archivos:
            self.archivos_sis = list(archivos)
            self.var_sis.set(f"{len(self.archivos_sis)} archivos cargados")

    def cargar_m1(self):
        archivos = filedialog.askopenfilenames(title="Molécula 1", filetypes=[("Archivos OUT", "*.out *.log")])
        if archivos:
            self.archivos_m1 = list(archivos)
            self.var_m1.set(f"{len(self.archivos_m1)} archivo(s) cargado(s)")

    def cargar_m2(self):
        archivos = filedialog.askopenfilenames(title="Molécula 2", filetypes=[("Archivos OUT", "*.out *.log")])
        if archivos:
            self.archivos_m2 = list(archivos)
            self.var_m2.set(f"{len(self.archivos_m2)} archivo(s) cargado(s)")

    def guardar_grafica(self):
        if not self.datos_grafica: return
        ruta = filedialog.asksaveasfilename(
            title="Guardar Gráfica",
            defaultextension=".png",
            filetypes=[("Imagen PNG", "*.png"), ("Todos los archivos", "*.*")]
        )
        if ruta:
            try:
                self.fig.savefig(ruta, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Guardado", f"Gráfica guardada exitosamente en:\n{ruta}")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo guardar la imagen:\n{e}")

    def procesar_y_graficar(self):
        if not self.archivos_sis or not self.archivos_m1 or not self.archivos_m2:
            messagebox.showwarning("Faltan Datos", "Debes cargar archivos para el Sistema, Molécula 1 y Molécula 2.")
            return
            
        try:
            atomos_m1 = int(self.spin_m1.get())
        except ValueError:
            return messagebox.showerror("Error", "Los átomos deben ser un número entero.")

        self.btn_calc.config(text="Procesando...", state="disabled")
        self.root.update()

        try:
            energias_m1 = [ExtractorGaussian.extraer_energia(f) for f in self.archivos_m1]
            energias_m2 = [ExtractorGaussian.extraer_energia(f) for f in self.archivos_m2]
            
            if None in energias_m1 or None in energias_m2:
                raise ValueError("Falta energía en archivos aislados.")

            e_m1_unica = energias_m1[0] if len(energias_m1) == 1 else None
            e_m2_unica = energias_m2[0] if len(energias_m2) == 1 else None

            self.datos_grafica = []
            HARTREE_TO_KCAL = 627.509

            for i, ruta_sis in enumerate(self.archivos_sis):
                e_sis = ExtractorGaussian.extraer_energia(ruta_sis)
                dist = ExtractorGaussian.extraer_distancia_centroides(ruta_sis, atomos_m1)
                
                if e_sis is None or dist is None: continue
                
                e_1 = e_m1_unica if e_m1_unica is not None else energias_m1[i]
                e_2 = e_m2_unica if e_m2_unica is not None else energias_m2[i]
                
                delta_e_kcal = (e_sis - (e_1 + e_2)) * HARTREE_TO_KCAL
                self.datos_grafica.append((dist, delta_e_kcal))

            if not self.datos_grafica:
                raise ValueError("No se extrajeron datos válidos.")

            # Ordenar por distancia
            self.datos_grafica.sort(key=lambda x: x[0])
            distancias = [d[0] for d in self.datos_grafica]
            energias = [d[1] for d in self.datos_grafica]

            # --- DIBUJAR GRÁFICA ---
            self.ax.clear()
            
            # 1. Trazar datos DFT (Gaussian)
            self.ax.plot(distancias, energias, marker='o', linestyle='-', color='#D32F2F', linewidth=2, markersize=6, label="DFT (Gaussian)")
            
            # 2. Trazar Lennard-Jones si está activo (Ajuste por Mínimos Cuadrados Robusto)
            if self.var_mostrar_lj.get():
                try:
                    # Definimos la función LJ con desplazamiento vertical (C)
                    def lennard_jones_desplazado(r, A, B, C):
                        return (A / r**12) - (B / r**6) + C
                    
                    r_data = np.array(distancias)
                    e_data = np.array(energias)
                    
                    # Limites: A y B deben ser positivos (0 a infinito), C puede ser cualquier valor
                    limites_inferiores = [0, 0, -np.inf]
                    limites_superiores = [np.inf, np.inf, np.inf]
                    
                    popt, _ = curve_fit(lennard_jones_desplazado, r_data, e_data, maxfev=100000, bounds=(limites_inferiores, limites_superiores))
                    A_opt, B_opt, C_opt = popt
                    
                    print(f"--- Parámetros LJ Encontrados ---")
                    print(f"A: {A_opt:.4e}")
                    print(f"B: {B_opt:.4e}")
                    print(f"C (Desplazamiento): {C_opt:.4f} kcal/mol")
                    
                    r_smooth = np.linspace(min(distancias) * 0.95, max(distancias) * 1.05, 200)
                    e_lj = lennard_jones_desplazado(r_smooth, A_opt, B_opt, C_opt)
                    
                    self.ax.plot(r_smooth, e_lj, linestyle='--', color='#1976D2', linewidth=2, label="Ajuste Teórico LJ")
                except Exception as e:
                    print(f"No se pudo realizar el ajuste de Lennard-Jones: {e}")

            # Decoración final
            self.ax.set_title("Curva de Energía de Interacción Molecular", fontsize=14, fontweight='bold', pad=10)
            self.ax.set_xlabel("Distancia entre Centroides (Å)", fontsize=11)
            self.ax.set_ylabel("Energía de Interacción ($\Delta E$) [kcal/mol]", fontsize=11)
            self.ax.grid(True, linestyle='--', alpha=0.7)
            self.ax.axhline(y=0, color='black', linestyle=':', alpha=0.6)
            self.ax.legend()
            
            # Ajustar el eje Y para mostrar absolutamente todos los puntos del sistema
            y_min = min(energias)
            y_max = max(energias)
            
            # Calculamos un margen del 5% para que los puntos no toquen los bordes de la ventana
            margen = (y_max - y_min) * 0.05 if y_max != y_min else 1.0
            self.ax.set_ylim(y_min - margen, y_max + margen)

            self.canvas.draw()
            
            # Habilitar el botón de guardado
            self.btn_guardar.config(state="normal")
            messagebox.showinfo("Éxito", f"Se graficaron {len(self.datos_grafica)} puntos.")

        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            self.btn_calc.config(text="CALCULAR Y GRAFICAR", state="normal")

if __name__ == "__main__":
    root = tk.Tk()
    app = AppGraficador(root)
    root.mainloop()
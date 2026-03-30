import os
import math
import tkinter as tk
from tkinter import filedialog, messagebox
# Importamos ttk para estilos modernos
from tkinter import ttk 

# ==========================================
# Núcleo Lógico (Backend)
# ==========================================
class SeparadorMoleculas:
    def __init__(self, archivo_gjf):
        self.archivo_gjf = archivo_gjf
        self.encabezado = []
        self.titulo = []
        self.carga_mult = ""
        self.coordenadas = []
        self.pie = []
        self._parsear_archivo()

    def _parsear_archivo(self):
        """Lee el archivo .gjf y separa sus secciones."""
        if not os.path.exists(self.archivo_gjf):
            return

        with open(self.archivo_gjf, 'r', encoding='utf-8') as f:
            lineas = f.readlines()

        seccion = 0
        for linea in lineas:
            linea_limpia = linea.strip()

            if linea_limpia == "" and seccion < 3:
                seccion += 1
                continue

            if seccion == 0:
                self.encabezado.append(linea)
            elif seccion == 1:
                self.titulo.append(linea)
            elif seccion == 2:
                if not self.carga_mult:
                    self.carga_mult = linea
                else:
                    partes = linea_limpia.split()
                    if len(partes) >= 4 and not partes[0].isdigit():
                        self.coordenadas.append({
                            'elemento': partes[0],
                            'x': float(partes[1]),
                            'y': float(partes[2]),
                            'z': float(partes[3])
                        })
                    else:
                        seccion = 3
                        self.pie.append(linea)
            elif seccion >= 3:
                self.pie.append(linea)

    def generar_separacion(self, atomos_mol1, pasos, distancia_paso, directorio_salida, nombre_base, progreso_callback=None):
        """Calcula centroides, vector unitario y genera n archivos."""
        if not os.path.exists(directorio_salida):
            os.makedirs(directorio_salida)

        total_atomos = len(self.coordenadas)
        if atomos_mol1 >= total_atomos or atomos_mol1 <= 0:
            raise ValueError(f"Número de átomos de Molécula 1 inválido (Total: {total_atomos})")

        # 1. Centroide Molécula 1
        mol1 = self.coordenadas[:atomos_mol1]
        c1 = [sum(a[eje] for a in mol1) / len(mol1) for eje in ['x', 'y', 'z']]

        # 2. Centroide Molécula 2
        mol2 = self.coordenadas[atomos_mol1:]
        c2 = [sum(a[eje] for a in mol2) / len(mol2) for eje in ['x', 'y', 'z']]

        # 3. Vector Director y Magnitud
        v = [c2[i] - c1[i] for i in range(3)]
        magnitud = math.sqrt(sum(coord**2 for coord in v))
        
        if magnitud == 0:
            raise ValueError("Los centroides de las moléculas coinciden. No se puede definir dirección.")

        # 4. Vector Unitario
        u = [coord / magnitud for coord in v]

        # Generación de archivos
        for paso in range(pasos + 1):
            nombre_archivo = f"{nombre_base}_{paso}.gjf"
            ruta_completa = os.path.join(directorio_salida, nombre_archivo)

            with open(ruta_completa, 'w', encoding='utf-8') as f:
                f.writelines(self.encabezado)
                f.write("\n")
                
                titulo_texto = self.titulo[0].strip() if self.titulo else "Generado"
                f.write(f"{titulo_texto} - Distancia +{distancia_paso * paso:.2f}A - Paso {paso}\n\n")
                
                f.write(self.carga_mult)
                
                for i, atomo in enumerate(self.coordenadas):
                    if i < atomos_mol1:
                        # M1 estática
                        x, y, z = atomo['x'], atomo['y'], atomo['z']
                    else:
                        # M2 se desplaza sobre el vector unitario
                        desplazamiento_total = distancia_paso * paso
                        x = atomo['x'] + (u[0] * desplazamiento_total)
                        y = atomo['y'] + (u[1] * desplazamiento_total)
                        z = atomo['z'] + (u[2] * desplazamiento_total)
                    
                    f.write(f" {atomo['elemento']:<2} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")
                
                f.write("\n")
                if self.pie:
                    f.writelines(self.pie)
                    f.write("\n")
            
            # Actualizar barra de progreso si se proporciona callback
            if progreso_callback:
                progreso_callback(paso, pasos)
                
        return True

# ==========================================
# Interfaz Gráfica (Frontend con ttk)
# ==========================================
class AppSeparador:
    def __init__(self, root):
        self.root = root
        self.root.title("Separador Geométrico Gaussian")
        self.root.geometry("620x520")
        
        # Configurar estilo moderno (Aplica tema del sistema operativo si es posible)
        self.style = ttk.Style()
        self.style.theme_use('clam') # 'clam', 'alt', 'default' son buenas opciones en Windows
        
        # Variables de control
        self.ruta_archivo = tk.StringVar()
        self.total_atomos = tk.IntVar(value=0)
        self.separador_logica = None

        self._crear_widgets()

    def _crear_widgets(self):
        # Margen general padding
        padd = {'padx': 15, 'pady': 8}
        
        # --- TITULO ---
        lbl_titulo = ttk.Label(self.root, text="Escaneo de Separación Molecular", font=("Helvetica", 16, "bold"))
        lbl_titulo.pack(pady=(20, 10))

        # --- SELECCIÓN DE ARCHIVO (FRAME) ---
        frame_archivo = ttk.LabelFrame(self.root, text=" 1. Archivo Plantilla .gjf ")
        frame_archivo.pack(fill="x", **padd)
        
        ttk.Entry(frame_archivo, textvariable=self.ruta_archivo, state="readonly").pack(side="left", fill="x", expand=True, padx=(10, 5), pady=10)
        ttk.Button(frame_archivo, text="Buscar...", command=self._seleccionar_archivo).pack(side="right", padx=(5, 10), pady=10)
        
        self.lbl_info_atomos = ttk.Label(frame_archivo, text="Átomos detectados: 0", font=("Helvetica", 9, "italic"))
        self.lbl_info_atomos.pack(pady=(0, 5), anchor="w", padx=15)

        # --- CONFIGURACIÓN DE MOLÉCULAS Y PASOS (FRAME) ---
        frame_config = ttk.LabelFrame(self.root, text=" 2. Parámetros de Simulación ")
        frame_config.pack(fill="x", **padd)
        
        # Grid para inputs
        grid_frame = ttk.Frame(frame_config)
        grid_frame.pack(pady=10, padx=10)
        
        # Inputs numéricos con Spinbox
        s_padd = {'padx': 10, 'pady': 5}
        
        ttk.Label(grid_frame, text="Átomos en Molécula 1:").grid(row=0, column=0, sticky="e", **s_padd)
        self.spin_m1 = ttk.Spinbox(grid_frame, from_=1, to=999, width=10)
        self.spin_m1.grid(row=0, column=1, sticky="w", **s_padd)
        
        ttk.Label(grid_frame, text="Cantidad de Pasos (n):").grid(row=1, column=0, sticky="e", **s_padd)
        self.spin_pasos = ttk.Spinbox(grid_frame, from_=1, to=500, width=10)
        self.spin_pasos.set(20) # Valor por defecto
        self.spin_pasos.grid(row=1, column=1, sticky="w", **s_padd)
        
        ttk.Label(grid_frame, text="Distancia por Paso (Å):").grid(row=1, column=2, sticky="e", **s_padd)
        self.entry_dist = ttk.Entry(grid_frame, width=12)
        self.entry_dist.insert(0, "0.25") # Valor por defecto
        self.entry_dist.grid(row=1, column=3, sticky="w", **s_padd)

        # --- SALIDA (FRAME) ---
        frame_salida = ttk.LabelFrame(self.root, text=" 3. Configuración de Salida ")
        frame_salida.pack(fill="x", **padd)
        
        ttk.Label(frame_salida, text="Nombre base archivos:").pack(side="left", padx=(15, 5), pady=10)
        self.entry_nombre = ttk.Entry(frame_salida)
        self.entry_nombre.insert(0, "separacion_scan")
        self.entry_nombre.pack(side="left", fill="x", expand=True, padx=(5, 15), pady=10)

        # --- BARRA DE PROGRESO Y BOTÓN EJECUTAR ---
        self.progress = ttk.Progressbar(self.root, orient="horizontal", mode="determinate")
        self.progress.pack(fill="x", padx=20, pady=(15, 5))
        
        self.btn_run = ttk.Button(self.root, text="GENERAR ARCHIVOS GJFs", command=self._ejecutar_proceso)
        # Estilo especial para el botón principal (si el tema lo soporta)
        self.style.configure('Action.TButton', font=('Helvetica', 11, 'bold'))
        self.btn_run.configure(style='Action.TButton')
        self.btn_run.pack(pady=15, ipady=5)

    def _seleccionar_archivo(self):
        file_path = filedialog.askopenfilename(
            title="Seleccionar plantilla Gaussian",
            filetypes=[("Archivos Gaussian", "*.gjf")]
        )
        if file_path:
            self.ruta_archivo.set(file_path)
            try:
                self.separador_logica = SeparadorMoleculas(file_path)
                n = len(self.separador_logica.coordenadas)
                self.total_atomos.set(n)
                self.lbl_info_atomos.config(text=f"Átomos detectados: {n}", foreground="green")
                # Actualizar el to_ del spinbox dinámicamente
                self.spin_m1.config(to=n-1) 
                if self.spin_m1.get() == "1" and n > 2: self.spin_m1.set(int(n/2))
            except Exception as e:
                messagebox.showerror("Error de Lectura", f"No se pudo parsear el archivo:\n{e}")
                self.ruta_archivo.set("")
                self.lbl_info_atomos.config(text="Error al leer archivo", foreground="red")

    def _actualizar_progreso(self, actual, total):
        porcentaje = (actual / total) * 100
        self.progress['value'] = porcentaje
        self.root.update_idletasks() # Forzar actualización visual

    def _ejecutar_proceso(self):
        # Validaciones iniciales
        if not self.separador_logica:
            messagebox.showwarning("Faltan datos", "Por favor, selecciona primero un archivo .gjf válido.")
            return
        
        try:
            m1 = int(self.spin_m1.get())
            pasos = int(self.spin_pasos.get())
            dist = float(self.entry_dist.get())
            nombre = self.entry_nombre.get().strip()
            
            if not nombre: raise ValueError("El nombre base no puede estar vacío.")
            if dist <= 0: raise ValueError("La distancia debe ser positiva.")

            # Preparar rutas
            origen = os.path.dirname(self.ruta_archivo.get())
            # Creamos carpeta única basada en el nombre para no mezclar escaneos
            carpeta_resultados = f"scan_{nombre}"
            salida = os.path.join(origen, carpeta_resultados)
            
            # Bloquear botón para evitar doble clic
            self.btn_run.config(state="disabled", text="Procesando...")
            self.progress['value'] = 0
            
            # Ejecutar lógica
            sep_exito = self.separador_logica.generar_separacion(
                m1, pasos, dist, salida, nombre, self._actualizar_progreso
            )
            
            if sep_exito:
                messagebox.showinfo("¡Éxito!", f"Se han generado {pasos+1} archivos correctamente.\n\nCarpeta: {carpeta_resultados}\nRuta: {origen}")

        except ValueError as ve:
            messagebox.showerror("Error de Entrada", f"Verifica los números ingresados:\n{ve}")
        except Exception as e:
            messagebox.showerror("Error Crítico", f"Ocurrió un error inesperado:\n{e}")
        finally:
            # Reincorporar estado original del botón
            self.btn_run.config(state="normal", text="GENERAR ARCHIVOS GJFs")
            self.progress['value'] = 0

if __name__ == "__main__":
    root = tk.Tk()
    app = AppSeparador(root)
    root.mainloop()
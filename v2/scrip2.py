import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import os

class GeneradorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Generador de Archivos Gaussian (Modificador de Z)")
        self.root.geometry("600x650")
        
        # --- Variables de Control ---
        self.ruta_archivo = tk.StringVar()
        self.texto_a_reemplazar = tk.StringVar(value="1.00000000")
        self.nombre_base_salida = tk.StringVar()
        
        self.z_inicial = tk.DoubleVar(value=1.0)
        self.z_final = tk.DoubleVar(value=6.0)
        self.paso = tk.DoubleVar(value=0.25)

        self._crear_interfaz()

    def _crear_interfaz(self):
        # Marco principal con márgenes
        main_frame = tk.Frame(self.root, padx=20, pady=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Título
        tk.Label(main_frame, text="Generador de Coordenadas Z", font=("Helvetica", 16, "bold")).pack(pady=(0, 20))

        # --- SECCIÓN 1: Selección de Archivo ---
        frame_archivo = tk.LabelFrame(main_frame, text="1. Archivo Plantilla (.gjf)", padx=10, pady=10)
        frame_archivo.pack(fill=tk.X, pady=5)

        entry_archivo = tk.Entry(frame_archivo, textvariable=self.ruta_archivo)
        entry_archivo.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        btn_buscar = tk.Button(frame_archivo, text="📂 Seleccionar Archivo", command=self.seleccionar_archivo)
        btn_buscar.pack(side=tk.RIGHT)

        # --- SECCIÓN 2: Configuración de Texto ---
        frame_config = tk.LabelFrame(main_frame, text="2. Configuración de Reemplazo", padx=10, pady=10)
        frame_config.pack(fill=tk.X, pady=10)

        # Grid para alinear etiquetas y entradas
        tk.Label(frame_config, text="Valor en plantilla a buscar:").grid(row=0, column=0, sticky="w", pady=5)
        tk.Entry(frame_config, textvariable=self.texto_a_reemplazar).grid(row=0, column=1, sticky="ew", padx=10)
        
        tk.Label(frame_config, text="Nombre base para salida:").grid(row=1, column=0, sticky="w", pady=5)
        tk.Entry(frame_config, textvariable=self.nombre_base_salida).grid(row=1, column=1, sticky="ew", padx=10)
        
        frame_config.columnconfigure(1, weight=1)

        # --- SECCIÓN 3: Parámetros Numéricos ---
        frame_nums = tk.LabelFrame(main_frame, text="3. Parámetros Numéricos", padx=10, pady=10)
        frame_nums.pack(fill=tk.X, pady=5)

        # Usamos grid para organizar Z inicial, final y salto
        tk.Label(frame_nums, text="Z Inicial:").grid(row=0, column=0, padx=5)
        tk.Entry(frame_nums, textvariable=self.z_inicial, width=10).grid(row=0, column=1, padx=5)

        tk.Label(frame_nums, text="Z Final:").grid(row=0, column=2, padx=5)
        tk.Entry(frame_nums, textvariable=self.z_final, width=10).grid(row=0, column=3, padx=5)

        tk.Label(frame_nums, text="Salto (Incremento):").grid(row=0, column=4, padx=5)
        tk.Entry(frame_nums, textvariable=self.paso, width=10).grid(row=0, column=5, padx=5)

        # --- BOTÓN DE ACCIÓN ---
        btn_generar = tk.Button(main_frame, text="GENERAR ARCHIVOS", bg="#4CAF50", fg="white", 
                                font=("Helvetica", 12, "bold"), height=2, command=self.procesar)
        btn_generar.pack(fill=tk.X, pady=20)

        # --- CONSOLA DE SALIDA ---
        tk.Label(main_frame, text="Registro de actividades:").pack(anchor="w")
        self.log_widget = scrolledtext.ScrolledText(main_frame, height=10, state='disabled', font=("Consolas", 9))
        self.log_widget.pack(fill=tk.BOTH, expand=True)

    def log(self, mensaje):
        """Escribe mensajes en el área de texto inferior"""
        self.log_widget.config(state='normal')
        self.log_widget.insert(tk.END, f"> {mensaje}\n")
        self.log_widget.see(tk.END)
        self.log_widget.config(state='disabled')
        self.root.update_idletasks() # Fuerza actualización visual

    def seleccionar_archivo(self):
        filename = filedialog.askopenfilename(
            title="Seleccionar archivo plantilla",
            filetypes=[("Archivos Gaussian", "*.gjf"), ("Archivos de texto", "*.txt"), ("Todos", "*.*")]
        )
        if filename:
            self.ruta_archivo.set(filename)
            # Sugerir un nombre base automáticamente basado en el archivo seleccionado
            base_name = os.path.splitext(os.path.basename(filename))[0]
            if not self.nombre_base_salida.get():
                self.nombre_base_salida.set(f"{base_name}_scan")
            self.log(f"Archivo seleccionado: {filename}")

    def procesar(self):
        # 1. Obtener datos y validar
        plantilla_path = self.ruta_archivo.get()
        target_val = self.texto_a_reemplazar.get()
        base_name = self.nombre_base_salida.get()

        if not plantilla_path or not os.path.exists(plantilla_path):
            messagebox.showerror("Error", "Por favor selecciona un archivo plantilla válido.")
            return
        
        if not base_name:
            messagebox.showerror("Error", "Debes especificar un nombre base para los archivos generados.")
            return

        try:
            inicio = self.z_inicial.get()
            final = self.z_final.get()
            salto = self.paso.get()
        except tk.TclError:
            messagebox.showerror("Error", "Asegúrate de que los valores numéricos (Z Inicial, Final, Salto) sean números válidos.")
            return

        # 2. Leer archivo
        self.log("--- Iniciando proceso ---")
        try:
            with open(plantilla_path, 'r', encoding='utf-8') as f:
                contenido_plantilla = f.read()
            
            # Advertencia si no encuentra el texto
            if target_val not in contenido_plantilla:
                resp = messagebox.askyesno("Advertencia", f"El texto '{target_val}' no se encontró en la plantilla.\n¿Deseas continuar generando archivos sin cambios?")
                if not resp:
                    self.log("Proceso cancelado por el usuario.")
                    return

            # 3. Bucle de generación
            directorio_salida = os.path.dirname(plantilla_path)
            valor_actual = inicio
            contador = 0
            
            # Usamos un pequeño margen de error para comparaciones de float
            while valor_actual <= final + (salto / 1000.0):
                nombre_archivo = f"{base_name}_{contador}.gjf"
                ruta_completa = os.path.join(directorio_salida, nombre_archivo)
                
                valor_fmt = f"{valor_actual:.8f}"
                
                # Reemplazo principal
                nuevo_contenido = contenido_plantilla.replace(target_val, valor_fmt)
                
                # Escribir archivo
                with open(ruta_completa, 'w', encoding='utf-8') as f_out:
                    f_out.write(nuevo_contenido)
                
                self.log(f"Generado: {nombre_archivo} [Z={valor_fmt}]")
                
                valor_actual += salto
                contador += 1
            
            self.log(f"--- ÉXITO: {contador} archivos creados en {directorio_salida} ---")
            messagebox.showinfo("Proceso Completado", f"Se han generado {contador} archivos correctamente.")

        except Exception as e:
            self.log(f"ERROR: {str(e)}")
            messagebox.showerror("Error Crítico", f"Ocurrió un error inesperado:\n{e}")

if __name__ == "__main__":
    # Configuración de la ventana raíz
    root = tk.Tk()
    app = GeneradorGUI(root)
    # Bucle principal de la interfaz
    root.mainloop()
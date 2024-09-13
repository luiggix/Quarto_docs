from PyPDF2 import PdfWriter
import glob

notebook_pdf = [(f) for f in glob.glob('*.pdf')]

print(notebook_pdf)

merger = PdfWriter()

for pdf in notebook_pdf:
    merger.append(pdf)

merger.write("intro_Python.pdf")
merger.close()
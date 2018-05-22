# ses_gen
### What is ses_gen
A program for triangulating analytically represented solvent-excluded surface. Developed as a bachelor thesis at *Faculty of informatics of Masaryk University in Brno, Czech Republic*.
### Execution
ses_gen can be used either as command line tool or be controlled by GUI.
The `classes/artifacts/ses_gen_jar/` folder contains `ses_gen.jar` executable. It has to be executed from the command line.
Command line options:

```
Usage: <main class> [options] Folder containing the SES analytic data to be 
      triangulated 
  Options:
    --edgelength, -e
      Average edge length of triangles in mesh
      Default: 0.3
    --help, -h
      Show this info
    -gui
      Use gui, view the resulting mesh directly.
      Default: false
    -obj
      Name of the OBJ text file the mesh should be exported into(optional)
    -stl
      Name of the STL text file the mesh should be exported into(optional)
    -v
      Print out various debug info
      Default: false
```
It is recommended to run the application with JVM configurated with 1GB of initial heap size, so that when triangulating big molecules the JVM does not have to do a lot of memory allocations.
### Example
`java -Xms1G -Xmx1G -jar ses_gen.jar -gui` executes ses_gen with JVM of 1G of initial and maximum heap size and allows the user to control the application with graphical user interface.

`java -Xms1G -Xmx1G -jar ses_gen.jar -e 0.4 -obj ./output.obj folder_with_molecule_data` generates triangular mesh for the molecule in `folder_with_molecule_data` with the average triangle edge length 0.4 and saves the mesh
in output.obj file and exits.

Example data can be found in `src/main/resources/proteins/` folder.

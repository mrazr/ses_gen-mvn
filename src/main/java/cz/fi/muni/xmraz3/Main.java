package cz.fi.muni.xmraz3;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.jogamp.common.os.Platform;
import cz.fi.muni.xmraz3.gui.MainPanel;
import cz.fi.muni.xmraz3.mesh.MeshGeneration;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class Main {
    public static void main(String[] args) {
        JCommander jc = JCommander.newBuilder().addObject(SesConfig.sesconfig).build();
        jc.parse(args);
        if (SesConfig.help){
            jc.usage();
            return;
        }
        if (SesConfig.inputFolder == null && (SesConfig.objFile != null || SesConfig.stlFile != null)){
            throw new ParameterException("Missing input folder argument");
        }
        if (SesConfig.objFile != null){
            Path p = Paths.get(SesConfig.objFile);
            if (Files.isDirectory(p)){
                SesConfig.objFile = p.resolve(Paths.get(SesConfig.inputFolder).getFileName().toString() + ".obj").toString();
            }
        }
        if (SesConfig.stlFile != null){
            Path p = Paths.get(SesConfig.stlFile);
            if (Files.isDirectory(p)){
                SesConfig.stlFile = p.resolve(Paths.get(SesConfig.inputFolder).getFileName().toString() + ".stl").toString();
            }
        }
        if (SesConfig.useGUI){
            MainPanel.startGUI(args);
        } else {
            if (SesConfig.inputFolder == null){
                jc.usage();
                return;
            }
            String raw = SurfaceParser.loadFile(Paths.get(SesConfig.inputFolder).resolve("info.json").toString());
            SurfaceParser.parseSesConfig(raw);
            SurfaceParser.ses_start(SesConfig.inputFolder);
            writeResults();
        }
    }

    private static void writeResults(){
        Runtime r = Runtime.getRuntime();
        long totalHeap = r.totalMemory();
        long usedHeap = totalHeap - r.freeMemory();
	float _usedHeap = (float)Math.ceil(usedHeap / (1024 * 1024));
        try (FileWriter fw = new FileWriter(Paths.get(SesConfig.inputFolder).getFileName().toString() + ".txt", true)){
            fw.write(Double.toString(SesConfig.edgeLimit));
            fw.write(" & ");
            fw.write(Long.toString(MeshGeneration.convexMeshTime));
            fw.write(" & ");
            fw.write(Long.toString(MeshGeneration.concaveMeshTime));
            fw.write(" & ");
            fw.write(Long.toString(MeshGeneration.toriMeshTime));
            fw.write(" & ");
            fw.write(Long.toString(MeshGeneration.trianglesGenerated.get()));
            fw.write(" & ");
            fw.write(Float.toString(_usedHeap));
            fw.write("\\\\");
            fw.write(Platform.getNewline());
        } catch (IOException e){
            e.printStackTrace();
        }
    }
}

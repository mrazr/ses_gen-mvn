package cz.fi.muni.xmraz3.mesh;

import cz.fi.muni.xmraz3.math.Point;
import cz.fi.muni.xmraz3.math.Vector;

import java.util.ArrayList;
import java.util.List;

public class ToroidalPatch {
    public Point probe1;
    public Point probe2;
    public Point midProbe;
    public List<Arc> convexPatchArcs;
    public List<Arc> concavePatchArcs;

    public int arcVertsCount;
    public CuspTriangle tr1;
    public CuspTriangle tr2;

    public List<Point> vertices;
    public List<Vector> normals;
    public int[] faces;
    public Point[] probes;
    public int vboOffset;
    public boolean circular = false;
    public boolean valid = true;
    public int id = -1;
    public static int nextID = 0;

    public ToroidalPatch(Point probe1, Point probe2, Point midProbe) {
        this.probe1 = probe1;
        this.probe2 = probe2;
        this.midProbe = midProbe;
        id = nextID++;

        convexPatchArcs = new ArrayList<>();
        concavePatchArcs = new ArrayList<>();
        vertices = new ArrayList<>();
        normals = new ArrayList<>();
    }
}

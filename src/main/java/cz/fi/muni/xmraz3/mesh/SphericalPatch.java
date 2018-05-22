package cz.fi.muni.xmraz3.mesh;

import cz.fi.muni.xmraz3.math.Point;
import cz.fi.muni.xmraz3.math.Sphere;
import cz.fi.muni.xmraz3.math.Vector;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class SphericalPatch {
    public Sphere sphere;
    public Vector patchNormal;
    public List<Arc> arcs;
    public List<Boundary> boundaries;
    public Map<Integer, List<ToroidalPatch>> tori;
    public List<SphericalPatch> neighbours;
    public List<Integer> intersectingPatches;
    public List<Point> vertices;
    public int[] faces;

    public int id;
    public int nextVertexID = 0;
    public static int nextConvexID = 0;
    public static int nextConcaveID = 0;
    public boolean convexPatch = true;
    public boolean valid = true;
    public boolean meshed = false;

    //opengl stuff
    public int lineCount = 0;
    public int vboOffset = 0;
    public int eboOffset = 0;
    public int lineOffset = 0;

    public SphericalPatch(Point center, double radius, boolean convex){
        this(new Sphere(center, radius), convex);
    }

    public SphericalPatch(Sphere s, boolean convex){
        this.sphere = s;
        id = (convex) ? nextConvexID++ : nextConcaveID++;
        boundaries = new ArrayList<>();
        arcs = new ArrayList<>();
        neighbours = new ArrayList<>();
        intersectingPatches = new ArrayList<>();
        tori = new TreeMap<>();
        vertices = new ArrayList<>();
        convexPatch = convex;
    }

    public SphericalPatch(Boundary b){
        boundaries = new ArrayList<>();
        boundaries.add(b);
    }
}

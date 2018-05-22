package cz.fi.muni.xmraz3.mesh;

import cz.fi.muni.xmraz3.SesConfig;
import cz.fi.muni.xmraz3.Surface;
import cz.fi.muni.xmraz3.math.Point;
import cz.fi.muni.xmraz3.math.Sphere;
import cz.fi.muni.xmraz3.math.Vector;
import cz.fi.muni.xmraz3.utils.ArcUtil;
import cz.fi.muni.xmraz3.utils.PatchUtil;

import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

public class MeshGeneration {
    private static int THREAD_COUNT = 0;
    public static MeshGeneration instance = new MeshGeneration();
    private static AtomicInteger threads_done = new AtomicInteger(0);

    public static AtomicInteger trianglesGenerated = new AtomicInteger(0);
    private static int[] _triangles;
    public static AtomicBoolean free = new AtomicBoolean(true);
    public static AtomicBoolean finished = new AtomicBoolean(false);
    private static Vector[] v;
    private static List<Face> facePool;
    private static AdvancingFrontMethod[] afms;

    public static long convexMeshTime = 0;
    public static long concaveMeshTime = 0;
    public static long toriMeshTime = 0;

    private MeshGeneration(){
        THREAD_COUNT = Runtime.getRuntime().availableProcessors();
        v = new Vector[THREAD_COUNT];

        _triangles = new int[THREAD_COUNT];
        for (int i = 0; i < THREAD_COUNT; ++i){
            _triangles[i] = 0;
            v[i] = new Vector(0, 0, 0);
        }
        facePool = new ArrayList<>(100);
        for (int i = 0; i < 100; ++i){
            facePool.add(i, new Face(0, 0, 0));
        }
    }

    private static void generateMesh(int start, int end, List<SphericalPatch> patches, int threadIdx) {
        free.set(false);
        AdvancingFrontMethod afm = afms[threadIdx];
        long startTime = System.currentTimeMillis();
        for (int i = start; i < end; ++i) {
            SphericalPatch a = patches.get(i);
            if (!a.valid){
                continue;
            }
            if (!a.meshed) {
                if (a.boundaries.size() > 0) {
                    ArcUtil.indexPoints(a);
                    afm.meshSphericalPatch(a);
                    if (afm.loop && SesConfig.verbose){
                        System.err.println((a.convexPatch) ? "convex " + i + " looped" : "concave " + i + " looped");
                    }
                    a.meshed = true;
                    afm.transferFacesToPatch();
                    _triangles[threadIdx] += a.faces.length / 3;
                }
            }
        }
        long endTime = System.currentTimeMillis();
        threads_done.incrementAndGet();
        if (threads_done.get() == THREAD_COUNT){
            if (SesConfig.verbose) {
                System.out.println(((patches.get(0).convexPatch) ? "Convex" : "Concave") + " patches meshed in " + (endTime - startTime) + " ms");
            }
            if (patches.get(0).convexPatch){
                convexMeshTime = (endTime - startTime);
            } else {
                concaveMeshTime = (endTime - startTime);
            }
            free.set(true);
            if (!patches.get(0).convexPatch){
                for (int i = 0; i < THREAD_COUNT; ++i){
                    trianglesGenerated.addAndGet(_triangles[i]);
                }
                System.out.println("Total number of triangles generated: " + trianglesGenerated.get());
                finished.set(true);
                for (int i = 0; i < afms.length; ++i){
                    afms[i] = null;
                }
            }
        }
    }

    public static void reset(){
        trianglesGenerated.set(0);
        finished.set(false);
        for (int i = 0; i < THREAD_COUNT; ++i){
            _triangles[i] = 0;
        }
       afms = new AdvancingFrontMethod[THREAD_COUNT];
        for (int i = 0; i < THREAD_COUNT; ++i){
            afms[i] = new AdvancingFrontMethod();
        }
    }

    public static void startMeshing(){
        finished.set(false);
        long startTime = System.currentTimeMillis();
        for (ToroidalPatch tp : Surface.rectangles){
            tp.vertices.clear();
            tp.normals.clear();
            meshToroidalPatch(tp);
            if (tp.faces != null) {
                _triangles[0] += tp.faces.length / 3;
            }
        }
        long endTime = System.currentTimeMillis();
        if (SesConfig.verbose) {
            System.out.println("Toroidal patches meshed in " + (endTime - startTime) + " milliseconds");
        }
        toriMeshTime = (endTime - startTime);
        MeshGeneration.threads_done.set(0);
        int step = SesConfig.atomCount / THREAD_COUNT;
        for (int i = 0; i < THREAD_COUNT; ++i){
            final int start = i;
            final int _step = step;
            Runnable r = new Runnable() {
                @Override
                public void run() {
                    MeshGeneration.generateMesh(start * _step, (start == THREAD_COUNT - 1) ? SesConfig.atomCount : (start + 1) * _step, Surface.convexPatches, start);
                }
            };
            (new Thread(r)).start();
        }
        while (!MeshGeneration.free.get()){}
        MeshGeneration.threads_done.set(0);
        step = SesConfig.trianglesCount / THREAD_COUNT;
        for (int i = 0; i < THREAD_COUNT; ++i){
            final int start = i;
            final int _step = step;
            Runnable r = new Runnable() {
                @Override
                public void run() {
                    MeshGeneration.generateMesh(start * _step, (start == THREAD_COUNT - 1) ? SesConfig.trianglesCount : (start + 1) * _step, Surface.triangles, start);
                }
            };
            (new Thread(r)).start();
        }
        while (!MeshGeneration.free.get()){}
    }

    private static List<Point> _top = new ArrayList<>(17);
    private static Arc _left = new Arc(new Point(0, 0, 0), 1.0);
    private static Arc _topL = new Arc(new Point(0, 0, 0), 1.0);
    private static Arc _right = new Arc(new Point(0, 0, 0), 1.0);

    private static Vector toProbe = new Vector(0, 0, 0);
    private static Vector atom1ToAtom2 = new Vector(0, 0, 0);
    private static Point centerOfRotation = new Point(0, 0, 0);
    private static Point centerOfTorus = new Point(0, 0, 0);
    private static Vector v1 = new Vector(0, 0, 0);
    private static Point currProbe = new Point(0, 0, 0);
    private static Point prevProbe = new Point(0, 0, 0);

    private static int nextFaceID = 0;
    private static List<Point> leftVArc = new ArrayList<>(17);
    private static List<Point> rightVArc = new ArrayList<>(17);
    public static void meshToroidalPatch(ToroidalPatch tp){
        try {
            if (!tp.circular){
                if (tp.concavePatchArcs.size() < 2 && (tp.tr1 == null || tp.tr2 == null)){
                    tp.valid = false;
                    return;
                } else if (tp.concavePatchArcs.size() == 2){
                    if (!tp.concavePatchArcs.get(0).valid || !tp.concavePatchArcs.get(1).valid){
                        tp.valid = false;
                        return;
                    }
                }
            }
            if (tp.tr1 != null){
                int subdLevel = (int)Math.max(ArcUtil.getSubdivisionLevel(tp.tr1.left), ArcUtil.getSubdivisionLevel(tp.tr1.right));
                _top.clear();
                for (int i = 0; i < tp.tr1.base.vrts.size(); ++i){
                    _top.add(tp.tr1.cuspPoint);
                }
                tp.probes = new Point[2 * tp.tr1.base.vrts.size()];
                int arcLev = ArcUtil.getSubdivisionLevel(tp.tr1.left);
                if (arcLev < subdLevel){
                    ArcUtil.refineArc(tp.tr1.left, SesConfig.edgeLimit, true, subdLevel - arcLev, false);
                }
                _left.center.change(tp.tr1.left.center);
                _left.radius = tp.tr1.left.radius;
                _left.vrts.clear();
                _left.vrts.addAll(tp.tr1.left.vrts);
                ArcUtil.reverseArc(_left, true);
                _topL.center.change(tp.tr1.cuspPoint);
                _topL.radius = 0;
                _topL.vrts.clear();
                _topL.vrts.addAll(_top);
                arcLev = ArcUtil.getSubdivisionLevel(tp.tr1.right);
                if (arcLev < subdLevel){
                   ArcUtil.refineArc(tp.tr1.right, SesConfig.edgeLimit, true, subdLevel - arcLev, false);
                }
                meshToroidalPatch(tp, tp.tr1.base, _topL, _left, tp.tr1.right, true);

                _top.clear();
                for (int i = 0; i < tp.tr2.base.vrts.size(); ++i){
                    _top.add(tp.tr2.cuspPoint);
                }
                subdLevel = Math.max(ArcUtil.getSubdivisionLevel(tp.tr2.left), ArcUtil.getSubdivisionLevel(tp.tr2.right));
                arcLev = ArcUtil.getSubdivisionLevel(tp.tr2.left);
                if (arcLev < subdLevel){
                    ArcUtil.refineArc(tp.tr2.left, SesConfig.edgeLimit, true, subdLevel - arcLev, false);
                }
                _left.center.change(tp.tr2.left.center);
                _left.radius = tp.tr2.left.radius;
                _left.vrts.clear();
                _left.vrts.addAll(tp.tr2.left.vrts);
                ArcUtil.reverseArc(_left, true);
                _topL.vrts.clear();
                _topL.vrts.addAll(_top);
                tp.arcVertsCount = tp.tr1.left.vrts.size();
                arcLev = ArcUtil.getSubdivisionLevel(tp.tr2.right);
                if (arcLev < subdLevel){
                    ArcUtil.refineArc(tp.tr2.right, SesConfig.edgeLimit, true, subdLevel - arcLev, false);
                }
                meshToroidalPatch(tp, tp.tr2.base, _topL, _left, tp.tr2.right, true);
                transferFacesToPatch(tp);
                Surface.toriFacesCount += tp.faces.length / 3;
                return;
            }
            Arc bottom = tp.convexPatchArcs.get(0);
            Arc left = null;
            Arc right = null;
            Arc top = tp.convexPatchArcs.get(1);

            atom1ToAtom2.changeVector(tp.convexPatchArcs.get(0).owner.sphere.center, tp.convexPatchArcs.get(1).owner.sphere.center).makeUnit();
            toProbe.changeVector(tp.probe1, tp.convexPatchArcs.get(0).owner.sphere.center);
            atom1ToAtom2.multiply(atom1ToAtom2.dotProduct(toProbe));
            centerOfRotation.assignTranslation(tp.convexPatchArcs.get(0).owner.sphere.center, atom1ToAtom2);
            if (tp.concavePatchArcs.size() == 0){
                if (Point.distance(tp.probe1, centerOfRotation) < SesConfig.probeRadius) {
                    centerOfTorus.assignTranslation(tp.probe1, v1.changeVector(tp.probe2, tp.probe1).multiply(0.5f));
                    double centerToCuspLength = Math.sqrt(Math.pow(SesConfig.probeRadius, 2) - Math.pow(Point.distance(tp.probe1, tp.probe2) / 2.f, 2));
                    Point bottomCusp = Point.translatePoint(centerOfTorus, v1.changeVector(bottom.owner.sphere.center, top.owner.sphere.center).makeUnit().multiply(centerToCuspLength));
                    Point topCusp = Point.translatePoint(centerOfTorus, v1.changeVector(top.owner.sphere.center, bottom.owner.sphere.center).makeUnit().multiply(centerToCuspLength));
                    Arc topForBottomRect = new Arc(bottom.center, bottom.radius);
                    for (int i = 0; i < bottom.vrts.size(); ++i) {
                        topForBottomRect.vrts.add(bottomCusp);
                    }
                    tp.probes = new Point[2 * bottom.vrts.size()];
                    Point newCenter = (Point.distance(bottom.end2, Sphere.getContactPoint(new Sphere(tp.probe1, SesConfig.probeRadius), bottom.owner.sphere)) < 0.0001) ? tp.probe1 : tp.probe2;
                    left = new Arc(newCenter, SesConfig.probeRadius);
                    left.vrts.add(bottom.end2);
                    Point mid;
                    Vector v;
                    left.vrts.add(bottomCusp);
                    left.setEndPoints(bottom.end2, bottomCusp, true);
                    ArcUtil.refineArc(left, SesConfig.edgeLimit, false,0, false);
                    int subdLevel = ArcUtil.getSubdivisionLevel(left);
                    newCenter = (Point.distance(left.center, tp.probe1) < 0.0001) ? tp.probe2 : tp.probe1;

                    right = new Arc(newCenter, SesConfig.probeRadius);
                    right.vrts.add(bottom.end1);
                    right.vrts.add(bottomCusp);
                    right.setEndPoints(bottom.end1, bottomCusp, true);
                    int numOfDivs = (int)(Math.log10(left.vrts.size() - 1) / Math.log10(2));
                    ArcUtil.refineArc(right, SesConfig.edgeLimit, true, subdLevel, false);
                    meshToroidalPatch(tp, bottom, topForBottomRect, left, right, false);

                    Arc bottomForTopRect = new Arc(top.center, top.radius);
                    for (int i = 0; i < top.vrts.size(); ++i){
                        bottomForTopRect.vrts.add(topCusp);
                    }
                    newCenter = (Point.distance(top.end2, Sphere.getContactPoint(new Sphere(tp.probe1, SesConfig.probeRadius), top.owner.sphere)) < 0.0001) ? tp.probe1 : tp.probe2;
                    left = new Arc(newCenter, SesConfig.probeRadius);
                    left.vrts.add(top.end2);
                    left.vrts.add(topCusp);
                    left.setEndPoints(top.end2, topCusp, true);
                    ArcUtil.refineArc(left, SesConfig.edgeLimit, true,subdLevel, false);
                    newCenter = (Point.subtractPoints(left.center, tp.probe1).sqrtMagnitude() < 0.0001) ? tp.probe2 : tp.probe1;
                    right = new Arc(newCenter, SesConfig.probeRadius);
                    right.vrts.add(top.end1);

                    right.vrts.add(topCusp);
                    right.setEndPoints(top.end1, topCusp, true);
                    numOfDivs = (int)(Math.log10(left.vrts.size() - 1) / Math.log10(2));
                    ArcUtil.refineArc(right, SesConfig.edgeLimit, true, subdLevel, false);
                    tp.arcVertsCount = left.vrts.size();
                    meshToroidalPatch(tp, top, bottomForTopRect, left, right, false);
                    transferFacesToPatch(tp);
                    Surface.toriFacesCount += tp.faces.length / 3;
                } else {
                    Point newCenter = (Point.distance(bottom.end2, Sphere.getContactPoint(new Sphere(tp.probe1, SesConfig.probeRadius), bottom.owner.sphere)) < 0.0001) ? tp.probe1 : tp.probe2;
                    _left.center.change(newCenter);
                    _left.radius = SesConfig.probeRadius;
                    _left.vrts.clear();
                    _left.vrts.add(bottom.end2);
                    _left.vrts.add(top.end1);
                    _left.setEndPoints(bottom.end2, top.end1, true);
                    ArcUtil.refineArc(_left, SesConfig.edgeLimit, false, 0, false);
                    newCenter = (Point.distance(_left.center, tp.probe1) < 0.0001) ? tp.probe2 : tp.probe1;

                    _right.center.change(newCenter);
                    _right.radius = SesConfig.probeRadius;
                    _right.vrts.clear();
                    _right.vrts.add(bottom.end1);
                    _right.vrts.add(top.end2);
                    _right.setEndPoints(bottom.end1, top.end2, true);
                    ArcUtil.refineArc(_right, SesConfig.edgeLimit, false, 0, false);
                    if (_right.vrts.size() != _left.vrts.size()){
                        System.out.println("weird");
                    }
                    tp.probes = new Point[bottom.vrts.size()];
                    tp.arcVertsCount = _left.vrts.size();
                    meshToroidalPatch(tp, bottom, top, _left, _right, false);
                    transferFacesToPatch(tp);
                    Surface.toriFacesCount += tp.faces.length / 3;
                }
            } else {
                toProbe.changeVector(tp.probe1, bottom.owner.sphere.center).makeUnit().multiply(bottom.owner.sphere.radius + SesConfig.probeRadius);
                atom1ToAtom2.changeVector(top.owner.sphere.center, bottom.owner.sphere.center).makeUnit();
                atom1ToAtom2.multiply(toProbe.dotProduct(atom1ToAtom2));
                double probeToRotationAx = -42;
                probeToRotationAx = PatchUtil.getProbeAxisDistance(tp.probe1, top.owner.sphere.center, bottom.owner.sphere.center);

                if (probeToRotationAx - SesConfig.probeRadius < 0.0){

                } else {
                    left = (Point.distance(bottom.end2, tp.concavePatchArcs.get(0).end2) < 0.0001) ? tp.concavePatchArcs.get(0) : tp.concavePatchArcs.get(1);
                    ArcUtil.reverseArc(left, true);
                    right = (left == tp.concavePatchArcs.get(0)) ? tp.concavePatchArcs.get(1) : tp.concavePatchArcs.get(0);
                    if (left.vrts.size() != right.vrts.size()){
                        Arc fewer = (left.vrts.size() > right.vrts.size()) ? right : left;
                        Arc more = (left == fewer) ? right : left;
                        int diff = ArcUtil.getSubdivisionLevel(more) - ArcUtil.getSubdivisionLevel(fewer);
                        ArcUtil.refineArc(fewer, SesConfig.edgeLimit, true, diff, false);
                    }
                    tp.probes = new Point[bottom.vrts.size()];
                    tp.arcVertsCount = left.vrts.size();
                    meshToroidalPatch(tp, bottom, top, left, right, false);
                    transferFacesToPatch(tp);
                    Surface.toriFacesCount += tp.faces.length;
                    ArcUtil.reverseArc(left, true);
                }
            }
        } catch (Exception e){
            System.out.println("for toroidal patch id: " + tp.id);
            tp.valid = false;
            e.printStackTrace();
        }
    }
    private static void meshToroidalPatch(ToroidalPatch tp, Arc bottom, Arc top, Arc left, Arc right, boolean special){
        try {
            int vertexOffset = tp.vertices.size();
            int arcLen = left.vrts.size();
            int probeOffset = (tp.vertices.size() > 0) ? bottom.vrts.size() : 0;
            leftVArc.clear();
            leftVArc.addAll(left.vrts);
            prevProbe.setAsMidpoint(left.center, left.center);
            tp.probes[probeOffset] = new Point(prevProbe);
            tp.vertices.addAll(leftVArc);
            for (int i = 1; i < bottom.vrts.size(); ++i) {
                Point vert = bottom.vrts.get(bottom.vrts.size() - i - 1);
                toProbe.changeVector(vert, bottom.owner.sphere.center).makeUnit().multiply(SesConfig.probeRadius);
                currProbe.assignTranslation(vert, toProbe);
                rightVArc.clear();
                if (i == bottom.vrts.size() - 1) {
                    rightVArc.addAll(right.vrts);

                    for (Point bod : rightVArc){
                        if (bod == null){
                            System.out.println(" ");
                        }
                    }
                } else {
                    if (!special) {
                        rightVArc.clear();
                        rightVArc = ArcUtil.generateCircArc(bottom.vrts.get(bottom.vrts.size() - i - 1), top.vrts.get(i), currProbe, SesConfig.probeRadius, left.vrts.size() - 1, false, rightVArc);
                        if (rightVArc.size() != left.vrts.size()){
                            System.out.println("incorrect number of vertices");
                        }
                    } else {
                        rightVArc.clear();
                        ArcUtil.generateCircArc(bottom.vrts.get(bottom.vrts.size() - i - 1), top.vrts.get(i), currProbe, SesConfig.probeRadius, left.vrts.size() - 1, false, rightVArc);
                    }
                }
                int l = (tp.vertices.size() - vertexOffset) / arcLen - 1;
                int r = l + 1;
                int m = arcLen;
                tp.vertices.addAll(rightVArc);
                tp.probes[r + probeOffset] = new Point(currProbe);
                for (int j = 0; j < left.vrts.size() - 1; ++j){
                    if (nextFaceID >= facePool.size()){
                         facePool.add(new Face(0, 0, 0));
                         facePool.add(new Face(0, 0, 0));
                     }
                     Face f1 = facePool.get(nextFaceID++);
                     Face f2 = facePool.get(nextFaceID++);
                     f1.a = l * m + j + vertexOffset;
                     f1.b = r * m + j + vertexOffset;
                     f1.c = l * m + j + vertexOffset + 1;

                     f2.a = l * m + j + vertexOffset + 1;
                     f2.b = r * m + j + vertexOffset;
                     f2.c = r * m + j + vertexOffset + 1;
                     Surface.numoftriangles += 2;
                }
                leftVArc.clear();
                leftVArc.addAll(rightVArc);
                prevProbe.setAsMidpoint(currProbe, currProbe);
            }

        } catch (Exception e){
            e.printStackTrace();
            tp.valid = false;
            System.err.println("Toroidal patch id: " + tp.id);
        }
    }

    private static void transferFacesToPatch(ToroidalPatch tp){
        tp.faces = new int[3 * nextFaceID];
        int j = 0;
        for (int i = 0; i < nextFaceID; ++i){
            Face f = facePool.get(i);
            tp.faces[j] = f.a;
            tp.faces[j + 1] = f.b;
            tp.faces[j + 2] = f.c;
            j += 3;
        }
        nextFaceID = 0;
    }

}

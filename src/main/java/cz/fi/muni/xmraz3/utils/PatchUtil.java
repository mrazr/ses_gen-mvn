package cz.fi.muni.xmraz3.utils;

import com.jogamp.opengl.math.Quaternion;
import cz.fi.muni.xmraz3.*;
import cz.fi.muni.xmraz3.math.Plane;
import cz.fi.muni.xmraz3.math.Point;
import cz.fi.muni.xmraz3.math.Sphere;
import cz.fi.muni.xmraz3.math.Vector;
import cz.fi.muni.xmraz3.mesh.*;
import smile.neighbor.Neighbor;

import java.util.*;
import java.util.function.Predicate;

public class PatchUtil {
    private static Vector v1 = new Vector(0, 0, 0);
    private static Vector v2 = new Vector(0, 0, 0);
    private static Vector v3 = new Vector(0, 0, 0);
    private static Plane p1 = new Plane(new Point(0, 0, 0), new Vector(0, 0, 0));
    private static Plane p2 = new Plane(new Point(0, 0, 0), new Vector(0, 0, 0));
    private static Vector circleN = new Vector(0, 0, 0);
    private static Point point = new Point(0, 0, 0);
    private static List<Point> middlevrts = new ArrayList<>(17);
    private static List<Arc> newArcs = new ArrayList<>();
    private static List<Boundary> newBs = new ArrayList<>();
    private static Plane _p = new Plane(new Point(0, 0, 0), new Vector(0, 0, 0));
    private static List<Boundary> newBS = new ArrayList<>();
    private static List<Boundary> toRemove = new ArrayList<>();
    private static List<Boundary> toRemove2 = new ArrayList<>();
    private static List<Plane> planePool = new ArrayList<>(50);
    private static boolean planePoolInitialized = false;
    private static int nextPlaneID = 0;
    private static Map<Integer, Map<Integer, List<Point>>> moip = new TreeMap<>();
    private static SphericalPatch curr;
    private static Plane currCirc;
    private static double currRad;
    private static List<Point> currInt;
    private static int currIter;
    private static Boundary _b = new Boundary();
    private static Arc _a1 = new Arc(new Point(0, 0, 0), 1.0);
    private static Arc _a2 = new Arc(new Point(0, 0, 0), 1.0);
    private static List<Point> vrtsPool = new ArrayList<>(17);
    private static List<Neighbor<double[], SphericalPatch>> neighbors = new ArrayList<>(50);
    private static Plane rho = new Plane(new Point(0, 0, 0), new Vector(0, 0, 0));
    private static List<Boundary> removeFromSP = new ArrayList<>();
    private static List<Boundary> processed = new ArrayList<>();
    private static Vector vInt = new Vector(0, 0, 0);
    private static Point pInt = new Point(0, 0, 0);
    private static Vector hypo = new Vector(0, 0, 0);
    private static Point in1 = new Point(0, 0, 0);
    private static Point in2 = new Point(0, 0, 0);
    private static Point midOfChord = new Point(0, 0, 0);
    private static Vector n = new Vector(0, 0, 0);
    private static List<Point> intersectionPoints = new ArrayList<>();
    private static List<Arc> exclude = new ArrayList<>();
    private static List<Point> invalid = new ArrayList<>();
    private static List<Point> usedPoints = new ArrayList<>();
    private static List<Point> usablePoints = new ArrayList<>();
    private static List<Point> usedPoints2 = new ArrayList<>();
    private static Point _nextPoint = new Point(0, 0, 0);
    private static Vector _genVector = new Vector(0, 0, 0);
    private static Quaternion _genQuaternion = new Quaternion();
    private static float[] _floatVector = new float[3];

    /*
    Method that treats a self-intersection of a given toroidal patch - trims its two spherical triangles against each other
    and assigns the toroidal patch with two CuspTriangle-s that are used later in the meshing phase to generate to two partial patches
     */
    public static void processSelfIntersectingTori(){
        for (ToroidalPatch tp : Surface.selfIntersectingRects){
            PatchUtil.torProcessSelfIntersection(tp);
        }
    }

    public static void torProcessSelfIntersection(ToroidalPatch tp){
        try {
            if (!tp.circular) {
                if (tp.concavePatchArcs.size() < 2) {
                    if (SesConfig.verbose) {
                        System.out.println("Incomplete rolling patch " + tp.id + ", skipping");
                    }
                    tp.valid = false;
                    return;
                }
            }
            if (!tp.concavePatchArcs.get(0).valid || !tp.concavePatchArcs.get(1).valid){
                tp.valid = false;
            }
            if (!tp.valid){
                Surface.rectangles.remove(tp);
                return;
            }
            Arc leftL = tp.concavePatchArcs.get(0);
            Arc rightL = tp.concavePatchArcs.get(1);
            Arc bottom = (Point.distance(leftL.end2, tp.convexPatchArcs.get(0).end2) < 0.001) ? tp.convexPatchArcs.get(0) : tp.convexPatchArcs.get(1);
            Arc top = (tp.convexPatchArcs.get(0) == bottom) ? tp.convexPatchArcs.get(1) : tp.convexPatchArcs.get(0);
            if (leftL.owner.intersectingPatches.contains(rightL.owner.id)){
                Arc cpl1 = null;
                Arc cpl2 = null;
                for (int i = 0; i < leftL.owner.boundaries.get(0).arcs.size(); ++i){
                    Arc a = leftL.owner.boundaries.get(0).arcs.get(i);
                    if (a.intersecting){
                        cpl1 = a;
                        break;
                    }
                }
                for (int i = 0; i < rightL.owner.boundaries.get(0).arcs.size(); ++i){
                    Arc a = rightL.owner.boundaries.get(0).arcs.get(i);
                    if (a.intersecting){
                        cpl2 = a;
                        break;
                    }
                }
                processIntersectingArcsOnPatch(cpl2);
                processIntersectingArcsOnPatch(cpl1);

                SphericalPatch left = cpl1.owner;
                SphericalPatch right = cpl2.owner;

                CuspTriangle tr1 = new CuspTriangle();
                CuspTriangle tr2 = new CuspTriangle();

                tr1.base = bottom;
                tr2.base = top;

                Arc target = null;
                for (int i = 0; i < left.boundaries.size(); ++i){
                    Boundary b = left.boundaries.get(i);
                    for (int j = 0; j < b.arcs.size(); ++j){
                        Arc a = b.arcs.get(j);
                        if (Point.distance(tr1.base.end2, a.end2) < 0.001){
                            target = a;
                            break;
                        }
                    }
                }
                tr1.left = target;
                tr1.cuspPoint = tr1.left.end1;
                target = null;
                for (int i = 0; i < right.boundaries.size(); ++i){
                    Boundary b = right.boundaries.get(i);
                    for (int j = 0; j < b.arcs.size(); ++j){
                        Arc a = b.arcs.get(j);
                        if (Point.distance(tr1.base.end1, a.end1) < 0.001){
                            target = a;
                            break;
                        }
                    }
                }
                tr1.right = target;
                target = null;
                for (int i = 0; i < right.boundaries.size(); ++i){
                    Boundary b = right.boundaries.get(i);
                    for (int j = 0; j < b.arcs.size(); ++j){
                        Arc a = b.arcs.get(j);
                        if (Point.distance(tr2.base.end2, a.end2) < 0.001){
                            target = a;
                            break;
                        }
                    }
                }
                tr2.left = target;
                tr2.cuspPoint = tr2.left.end1;
                target = null;
                for (int i = 0; i < left.boundaries.size(); ++i){
                    Boundary b = left.boundaries.get(i);
                    for (int j = 0; j < b.arcs.size(); ++j){
                        Arc a = b.arcs.get(j);
                        if (Point.distance(tr2.base.end1, a.end1) < 0.001){
                            target = a;
                            break;
                        }
                    }
                }
                tr2.right = target;
                tp.tr1 = tr1;
                tp.tr2 = tr2;
                return;
            }
            Point circle = new Point(0, 0, 0);
            double radius = computeIntersectionCircle(leftL.owner.sphere.center, rightL.owner.sphere.center, circle, SesConfig.probeRadius);
            circleN.changeVector(leftL.owner.sphere.center, circle).makeUnit();
            Point[] cusps = new Point[2];
            cusps[0] = computeCusp(tp.probe1, tp.convexPatchArcs.get(0).owner.sphere, tp.convexPatchArcs.get(1).owner.sphere);
            cusps[1] = computeCusp(tp.probe1, tp.convexPatchArcs.get(1).owner.sphere, tp.convexPatchArcs.get(0).owner.sphere);

            if (cusps[0] == null || cusps[1] == null){
                return;
            }
            if (Point.distance(leftL.end1, cusps[0]) - Point.distance(leftL.end1, cusps[1]) > 0.0) {
                Point tmp = cusps[1];
                cusps[1] = cusps[0];
                cusps[0] = tmp;
            }
            Arc leftStart = new Arc(leftL.center, leftL.radius);
            leftStart.owner = leftL.owner;
            leftStart.setEndPoints(leftL.end1, cusps[0], true);
            leftStart.prev = leftL.prev;
            leftStart.prev.next = leftStart;
            leftStart.vrts.add(leftStart.end1);
            leftStart.vrts.add(cusps[0]);
            ArcUtil.refineArc(leftStart, Surface.maxEdgeLen, false, 0, false);

            Arc leftEnd = new Arc(leftL.center, leftL.radius);
            leftEnd.owner = leftL.owner;
            leftEnd.setEndPoints(cusps[1], leftL.end2, true);
            leftEnd.next = leftL.next;
            leftEnd.next.prev = leftEnd;
            leftEnd.vrts.add(cusps[1]);
            leftEnd.vrts.add(leftEnd.end2);
            ArcUtil.refineArc(leftEnd, Surface.maxEdgeLen, false, 0, false);

            Arc rightStart = new Arc(rightL.center, rightL.radius);
            rightStart.owner = rightL.owner;
            rightStart.setEndPoints(rightL.end1, cusps[1], true);
            rightStart.prev = rightL.prev;
            rightStart.prev.next = rightStart;
            rightStart.vrts.add(rightStart.end1);
            rightStart.vrts.add(cusps[1]);
            ArcUtil.refineArc(rightStart, Surface.maxEdgeLen, false,0, false);

            Arc rightEnd = new Arc(rightL.center, rightL.radius);
            rightEnd.owner = rightL.owner;
            rightEnd.setEndPoints(cusps[0], rightL.end2, true);
            rightEnd.next = rightL.next;
            rightEnd.next.prev = rightEnd;
            rightEnd.vrts.add(cusps[0]);
            rightEnd.vrts.add(rightEnd.end2);
            ArcUtil.refineArc(rightEnd, Surface.maxEdgeLen, false,0, false);

            v1.changeVector(cusps[0], circle).makeUnit(); //toS
            v2.changeVector(cusps[1], circle).makeUnit(); //toE
            n.assignNormalVectorOf(v2, v1).makeUnit();
            middlevrts.clear();
            if (n.dotProduct(circleN) > 0.0){
                ArcUtil.generateCircArc(cusps[0], cusps[1], circle, radius, 2, true, middlevrts);
            } else {
                ArcUtil.generateCircArc(cusps[0], cusps[1], circle, radius, 2, false, middlevrts);
            }
            Arc leftMid = new Arc(circle, radius);
            leftMid.owner = leftL.owner;
            leftMid.setEndPoints(cusps[0], cusps[1], false);
            leftMid.setNormal(v1.changeVector(leftL.owner.sphere.center, circle).makeUnit());
            leftMid.prev = leftStart;
            leftStart.next = leftMid;
            leftMid.next = leftEnd;
            leftEnd.prev = leftMid;
            leftMid.vrts.addAll(middlevrts);
            ArcUtil.refineArc(leftMid, Surface.maxEdgeLen, false,0, false);

            Arc rightMid = ArcUtil.cloneArc(leftMid);
            rightMid.owner = rightL.owner;
            rightMid.setNormal(leftMid.normal);
            ArcUtil.reverseArc(rightMid, true);
            rightMid.baseSubdivision = leftMid.baseSubdivision;
            rightMid.next = rightEnd;
            rightEnd.prev = rightMid;
            rightMid.prev = rightStart;
            rightStart.next = rightMid;
            rightMid.opposite = leftMid;
            leftMid.opposite = rightMid;

            rightEnd.vrts.remove(0);
            rightEnd.vrts.add(0, rightMid.end2);
            rightEnd.end1 = rightMid.end2;

            rightStart.vrts.remove(rightStart.vrts.size() - 1);
            rightStart.vrts.add(rightMid.end1);
            rightStart.end2 = rightMid.end1;

            SphericalPatch left = leftL.owner;
            SphericalPatch right = rightL.owner;
            /*if (left.intersectingPatches.contains(right.id) || right.intersectingPatches.contains(left.id)){
                System.out.println("found double trimming " + left.id + " " + right.id);
            }*/

            left.neighbours.add(right);
            right.neighbours.add(left);

            leftStart.owner = leftMid.owner = leftEnd.owner = left;
            rightStart.owner = rightMid.owner = rightEnd.owner = right;
            leftMid.intersecting = true;
            rightMid.intersecting = true;
            Surface.intersectingArcs.add(leftMid);
            Surface.intersectingArcs.add(rightMid);
            left.intersectingPatches.add(right.id);
            right.intersectingPatches.add(left.id);

            left.boundaries.get(0).arcs.clear();
            leftStart.bOwner = leftMid.bOwner = leftEnd.bOwner = left.boundaries.get(0);
            rightStart.bOwner = rightEnd.bOwner = rightMid.bOwner = right.boundaries.get(0);

            Arc l = leftStart;
            do {
                left.boundaries.get(0).arcs.add(l);
                l = l.next;
            } while (l != leftStart);
            right.boundaries.get(0).arcs.clear();
            l = rightStart;
            do {
                right.boundaries.get(0).arcs.add(l);
                l = l.next;
            } while (l != rightStart);
            left.boundaries.get(0).vrts.clear();
            for (Arc lo : left.boundaries.get(0).arcs){
                for (int i = 0; i < lo.vrts.size() - 1; ++i){
                    left.boundaries.get(0).vrts.add(lo.vrts.get(i));
                }
            }
            right.boundaries.get(0).vrts.clear();
            for (Arc lo : right.boundaries.get(0).arcs){
                for (int i = 0; i < lo.vrts.size() - 1; ++i){
                    right.boundaries.get(0).vrts.add(lo.vrts.get(i));
                }
            }
            ArcUtil.buildEdges(left.boundaries.get(0), false);
            ArcUtil.buildEdges(right.boundaries.get(0), false);

            tp.tr1 = new CuspTriangle();
            tp.tr2 = new CuspTriangle();
            tp.tr1.base = (Point.distance(leftStart.end1, tp.convexPatchArcs.get(0).end1) < 0.0001) ? tp.convexPatchArcs.get(0) : tp.convexPatchArcs.get(1);
            tp.tr1.left = rightEnd;
            tp.tr1.right = leftStart;
            tp.tr1.left.opposite = tp.tr1.right;
            tp.tr1.right.opposite = tp.tr1.left;
            rightEnd.cuspTriangle = leftStart.cuspTriangle = tp.tr1;
            if (rightEnd.vrts.size() != leftStart.vrts.size()){
                System.out.println(" ");
            }
            tp.tr1.cuspPoint = cusps[0];

            tp.tr2.base = (tp.tr1.base == tp.convexPatchArcs.get(0)) ? tp.convexPatchArcs.get(1) : tp.convexPatchArcs.get(0);
            tp.tr2.left = leftEnd;
            tp.tr2.right = rightStart;
            tp.tr2.left.opposite = tp.tr2.right;
            tp.tr2.right.opposite = tp.tr2.left;
            leftEnd.cuspTriangle = rightStart.cuspTriangle = tp.tr2;
            if (leftEnd.vrts.size() != rightStart.vrts.size()){
                System.out.println(" ");
            }
            tp.tr2.cuspPoint = cusps[1];
            leftStart.torus = leftEnd.torus = rightStart.torus = rightEnd.torus = null;
            tp.concavePatchArcs.clear();
        } catch (Exception e){
            System.err.println("tp id: " + tp.id);
            e.printStackTrace();
        }
    }
    /*
    After the treatment of self-intersecting toroidal patches, some spherical triangles might have arcs that intersect other arcs-
    they trimmed here.
     */
    public static void processSelfIntersectingConcavePatches() {
        for (Arc cpl : Surface.intersectingArcs) {
            PatchUtil.trimSelfIntersectingPatch(cpl);
        }
    }

    public static void trimSelfIntersectingPatch(Arc arc){
        try {
            if (!arc.valid || !arc.intersecting){
                return;
            }
            SphericalPatch sp = arc.owner;
            intersectionPoints.clear();
            exclude.clear();
            exclude.add(arc);
            invalid.clear();
            findIntersectionPoints(sp, arc.center, arc.radius, intersectionPoints, exclude);
            for (int i = 0; i < intersectionPoints.size(); ++i){
                Point v = intersectionPoints.get(i);
                if (Point.distance(v, arc.prev.end1) < 0.0015 || Point.distance(v, arc.next.end2) < 0.0015){
                    invalid.add(v);
                }
            }
            intersectionPoints.removeAll(invalid);
            if (intersectionPoints.size() > 1){ //theoretically at most 2 intersection points should be encountered(and in most cases), 1 point is special case which is not handled right now(it is rare)
                for (int i = 0; i < sp.boundaries.size(); ++i){
                    Boundary b = sp.boundaries.get(i);
                    for (int j = 0; j < b.arcs.size(); ++j){
                        Arc a = b.arcs.get(j);
                        a.valid = false;
                    }
                }
                Plane circle = rho;// new Plane(arc.center, arc.normal);
                rho.redefine(arc.center, arc.normal);
                newBS.clear();
                Point first = ArcUtil.findClosestPointOnCircle(intersectionPoints, arc.end1, true, arc.center, arc.normal, true);
                intersectionPoints.remove(first);
                Arc a = first.arc;//ArcUtil.findContainingArc(first, circle, sp, arc);
                if (Point.distance(a.end1, first) < 0.001 && nextSign(first, a, circle) < 0.0){
                    a = a.prev;
                } else if (Point.distance(a.end2, first) < 0.001 && nextSign(first, a, circle) > 0.0){
                    a = a.next;
                }
                if (a == null){
                    int c = 4;
                }
                Boundary b = new Boundary();
                Arc newA = (Point.distance(first, arc.end1) < 0.003) ? null : new Arc(arc.center, arc.radius);
                if (newA != null) {
                    newA.vrts.add(arc.end1);
                    newA.vrts.add(first);
                    newA.setEndPoints(arc.end1, first, false);
                    newA.setNormal(arc.normal);
                    if (arc.cuspTriangle != null){
                        int subdLevel = ArcUtil.getSubdivisionLevel(arc);
                        newA.cuspTriangle = arc.cuspTriangle;
                        if (newA.cuspTriangle.left == arc){
                            newA.cuspTriangle.left = newA;
                        } else {
                            newA.cuspTriangle.right = newA;
                        }
                        ArcUtil.refineArc(newA, 0, true, subdLevel, false);
                    } else {
                        ArcUtil.refineArc(newA, Surface.maxEdgeLen, false, 0, false);
                    }
                    if (arc.torus != null){
                        if (arc.torus.tr1 != null){
                            if (arc == arc.torus.tr1.left){
                                arc.torus.tr1.left = newA;
                            } else if (arc == arc.torus.tr1.right){
                                arc.torus.tr1.right = newA;
                            } else if (arc == arc.torus.tr2.left){
                                arc.torus.tr2.left = newA;
                            } else if (arc == arc.torus.tr2.right){
                                arc.torus.tr2.right = newA;
                            }
                        } else {
                            if (arc.torus.concavePatchArcs.get(0) == arc){
                                arc.torus.concavePatchArcs.remove(arc);
                                arc.torus.concavePatchArcs.add(newA);
                            } else {
                                arc.torus.concavePatchArcs.remove(arc);
                                arc.torus.concavePatchArcs.add(newA);
                            }
                        }
                        newA.torus = arc.torus;
                    }
                } else {
                    first = arc.end1;
                }
                Arc newA2 = new Arc(a.center, a.radius);
                newA2.vrts.add(first);
                newA2.vrts.add(a.end2);
                newA2.setEndPoints(first, a.end2, false);
                newA2.setNormal(a.normal);
                if (a.cuspTriangle != null){
                    int subdLevel = ArcUtil.getSubdivisionLevel(a);
                    newA2.cuspTriangle = a.cuspTriangle;
                    if (a.cuspTriangle.left == a){
                        newA2.cuspTriangle.left = newA2;
                    } else {
                        newA2.cuspTriangle.right = newA2;
                    }
                    ArcUtil.refineArc(newA2, 0, true, subdLevel, false);
                } else {
                    ArcUtil.refineArc(newA2, Surface.maxEdgeLen, false, 0, false);
                }
                if (a.torus != null){
                    if (a.torus.tr1 != null){
                        if (a == a.torus.tr1.left){
                            a.torus.tr1.left = newA2;
                        } else if (a == a.torus.tr1.right){
                            a.torus.tr1.right = newA2;
                        } else if (a == a.torus.tr2.left){
                            a.torus.tr2.left = newA2;
                        } else if (a == a.torus.tr2.right){
                            a.torus.tr2.right = newA2;
                        }
                    } else {
                        if (a.torus.concavePatchArcs.get(0) == a){
                            a.torus.concavePatchArcs.remove(a);
                            a.torus.concavePatchArcs.add(newA2);
                        } else {
                            a.torus.concavePatchArcs.remove(a);
                            a.torus.concavePatchArcs.add(newA2);
                        }
                    }
                    newA2.torus = a.torus;
                }
                if (newA != null) {
                    b.arcs.add(newA);
                }
                b.arcs.add(newA2);
                Arc start = arc;
                Arc a_ = a.next;
                while (a_ != start) {
                    b.arcs.add(a_);
                    if (a_ == null){
                        System.out.println(" ");
                    }
                    if (a_.next == null){
                        System.out.println(" ");
                    }
                    a_ = a_.next;
                }

                b.patch = sp;
                ArcUtil.buildEdges(b, true);

                newBS.add(b);

                Point second = intersectionPoints.get(0); //remaining intersection point
                a = second.arc;//ArcUtil.findContainingArc(second, circle, sp, arc);
                if (Point.distance(second, a.end1) < 0.001 && nextSign(second, a, circle) < 0.0){
                    a = a.prev;
                } else if (Point.distance(second, a.end2) < 0.001 && nextSign(second, a, circle) > 0.0){
                    a = a.next;
                }
                b = new Boundary();
                newA = new Arc(a.center, a.radius);
                newA.vrts.add(a.end1);
                newA.vrts.add(second);
                newA.setEndPoints(a.end1, second, false);
                newA.setNormal(a.normal);
                if (a.cuspTriangle != null){
                    int subdLevel = ArcUtil.getSubdivisionLevel(a);
                    newA.cuspTriangle = a.cuspTriangle;
                    if (a.cuspTriangle.left == a){
                        newA.cuspTriangle.left = newA;
                    } else {
                        newA.cuspTriangle.right = newA;
                    }
                    ArcUtil.refineArc(newA, 0, true, subdLevel, false);
                } else {
                    ArcUtil.refineArc(newA, Surface.maxEdgeLen, false, 0, false);
                }

                if (a.torus != null){
                    if (a.torus.tr1 != null){
                        if (a == a.torus.tr1.left){
                            a.torus.tr1.left = newA;
                        } else if (a == a.torus.tr1.right){
                            a.torus.tr1.right = newA;
                        } else if (a == a.torus.tr2.left){
                            a.torus.tr2.left = newA;
                        } else if (a == a.torus.tr2.right){
                            a.torus.tr2.right = newA;
                        }
                    } else {
                        if (a.torus.concavePatchArcs.get(0) == a){
                            a.torus.concavePatchArcs.remove(a);
                            a.torus.concavePatchArcs.add(newA);
                        } else {
                            a.torus.concavePatchArcs.remove(a);
                            a.torus.concavePatchArcs.add(newA);
                        }
                    }
                    newA.torus = a.torus;
                }
                newA2 = (Point.distance(second, arc.end2) < 0.003) ? null : new Arc(arc.center, arc.radius);
                if (newA2 != null) {
                    newA2.vrts.add(second);
                    newA2.vrts.add(arc.end2);
                    newA2.setEndPoints(second, arc.end2, false);
                    newA2.setNormal(arc.normal);
                    if (arc.cuspTriangle != null){
                        int subdLevel = ArcUtil.getSubdivisionLevel(arc);
                        newA2.cuspTriangle = arc.cuspTriangle;
                        if (arc.cuspTriangle.left == arc){
                            newA2.cuspTriangle.left = newA2;
                        } else {
                            newA2.cuspTriangle.right = newA2;
                        }
                        ArcUtil.refineArc(newA2, 0, true, subdLevel, false);
                    } else {
                        ArcUtil.refineArc(newA2, Surface.maxEdgeLen, false, 0, false);
                    }
                    if (arc.torus != null){
                        if (arc.torus.tr1 != null){
                            if (arc == arc.torus.tr1.left){
                                arc.torus.tr1.left = newA2;
                            } else if (arc == arc.torus.tr1.right){
                                arc.torus.tr1.right = newA2;
                            } else if (arc == arc.torus.tr2.left){
                                arc.torus.tr2.left = newA2;
                            } else if (arc == arc.torus.tr2.right){
                                arc.torus.tr2.right = newA2;
                            }
                        } else {
                            if (arc.torus.concavePatchArcs.get(0) == arc){
                                arc.torus.concavePatchArcs.remove(arc);
                                arc.torus.concavePatchArcs.add(newA2);
                            } else {
                                arc.torus.concavePatchArcs.remove(arc);
                                arc.torus.concavePatchArcs.add(newA2);
                            }
                        }
                        newA2.torus = arc.torus;
                    }
                } else {
                    second = arc.end2;
                    newA.end2 = second;
                    newA.vrts.remove(newA.vrts.size() - 1);
                    newA.vrts.add(second);
                }
                start = a;
                a_ = arc.next;
                b.arcs.add(newA);
                if (newA2 != null) {
                    b.arcs.add(newA2);
                }
                while (a_ != start) {
                    b.arcs.add(a_);
                    a_ = a_.next;
                }
                b.patch = sp;
                ArcUtil.buildEdges(b, true);
                newBS.add(b);
                sp.boundaries.clear();
                sp.boundaries.addAll(newBS);
            }
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void processIntersectingArcsOnPatch(Arc arc){
        if (!arc.intersecting || !arc.valid){
            return;
        }
        intersectionPoints.clear();
        Point lastPoI = null;
        boolean found = false;
        SphericalPatch sp = arc.owner;
        Vector u1 = v1.changeVector(arc.end1, arc.center).makeUnit();
        Vector u2 = v2.changeVector(arc.end2, arc.center).makeUnit();
        Vector n1 = v3.changeVector(sp.sphere.center, arc.center).makeUnit();
        p1.redefine(arc.center, n1);
        for (int i = 0; i < sp.boundaries.size(); ++i){
            Boundary b = sp.boundaries.get(i);
            for (int j = 0; j < b.arcs.size(); ++j){
                Arc k = b.arcs.get(j);
                if (k == arc || k == arc.prev || k == arc.next){
                    continue;
                }
                boolean allInside = true;
                for (int l = 0; l < k.vrts.size(); ++l){
                    if (!(p1.checkPointLocation(k.vrts.get(l)) > 0.0)){
                        allInside = false;
                        break;
                    }
                }
                if (allInside){
                    continue;
                }
                v1.changeVector(k.end1, k.center).makeUnit();
                v2.changeVector(k.end2, k.center).makeUnit();
                p2.redefine(k.center, v1, v2);
                if (Plane.getIntersectionLine(p1, p2, v1, pInt)){
                    hypo.changeVector(arc.center, pInt);
                    v1.makeUnit();
                    double odvesna = v1.dotProduct(hypo);
                    double dist = Math.sqrt(hypo.dotProduct(hypo) - Math.pow(odvesna, 2));
                    if (dist - arc.radius < 0.0){
                        midOfChord.assignTranslation(pInt, v1.multiply(odvesna));
                        double odv2 = Math.sqrt(Math.pow(arc.radius, 2) - Math.pow(dist, 2));
                        in1.assignTranslation(midOfChord, v1.makeUnit().multiply(odv2));
                        in2.assignTranslation(midOfChord, v1.multiply(-1.0));
                        if ((k.isInside(in1) && arc.isInside(in1)) || (k.isInside(in2) && arc.isInside(in2))){
                            found = true;
                        }
                        if (k.isInside(in1) && arc.isInside(in1)){
                            boolean foundSimilar = false;
                            for (int l = 0; l < intersectionPoints.size(); ++l){
                                if (Point.distance(in1, intersectionPoints.get(l)) < 0.001){
                                    foundSimilar = true;
                                    break;
                                }
                            }
                            if (!foundSimilar){
                                intersectionPoints.add(new Point(in1));
                            }
                        }
                        if (k.isInside(in2) && arc.isInside(in2)){
                            boolean foundSimilar = false;
                            for (int l = 0; l < intersectionPoints.size(); ++l){
                                if (Point.distance(in2, intersectionPoints.get(l)) < 0.001){
                                    foundSimilar = true;
                                    break;
                                }
                            }
                            if (!foundSimilar){
                                intersectionPoints.add(new Point(in2));
                            }
                        }
                    }
                }
            }
        }
        if (found && intersectionPoints.size() == 2) {
            for (int i = 0; i < sp.boundaries.size(); ++i){
                Boundary b = sp.boundaries.get(i);
                for (int j = 0; j < b.arcs.size(); ++j){
                    Arc a = b.arcs.get(j);
                    a.valid = false;
                }
            }
            newBs.clear();
            Boundary b = new Boundary();
            b.patch = sp;
            Arc start = arc.next.next;
            Arc a = start;
            do {
                _p.redefine(a.center, a.normal);
                boolean trim = false;
                Point in = null;
                for (int j = 0; j < intersectionPoints.size(); ++j){
                    Point i = intersectionPoints.get(j);
                    if (Math.abs(_p.checkPointLocation(i)) < 0.001 && a.isInside(i)){
                        trim = true;
                        if (in == null || Point.distance(in, a.end1) - Point.distance(i, a.end1) > 0.0){
                            in = i;
                        }
                    }
                }
                if (trim){
                    Arc newA2 = new Arc(a.center, a.radius);
                    newA2.setEndPoints(a.end1, in, false);
                    newA2.setNormal(a.normal);
                    newA2.vrts.add(a.end1);
                    newA2.vrts.add(in);
                    newA2.valid = true;
                    ArcUtil.refineArc(newA2, Surface.maxEdgeLen, false,0, false);
                    Arc newA = new Arc(arc.center, arc.radius);
                    newA.setNormal(arc.normal);
                    newA.vrts.add(in);
                    newA.vrts.add(arc.end2);
                    newA.setEndPoints(in, arc.end2, false);
                    ArcUtil.refineArc(newA, Surface.maxEdgeLen, false,0, false);

                    newA.prev = newA2;
                    newA2.next = newA;
                    newA.next = arc.next;
                    newA.next.prev = newA;
                    b.arcs.add(newA2);
                    b.arcs.add(newA);
                    a.intersecting = false;
                    a = newA;
                    a.intersecting = false;
                    a.valid = true;
                    intersectionPoints.remove(in);
                } else {
                    a.valid = true;
                    b.arcs.add(a);
                }
                a = a.next;
            } while (a != start);
            ArcUtil.buildEdges(b, true);
            newBs.add(b);
            b = new Boundary();
            b.patch = sp;
            newArcs.clear();
            start = arc.prev.prev;
            a = start;
            do {
                _p.redefine(a.center, a.normal);
                boolean trim = false;
                Point in = null;
                for (int j = 0; j < intersectionPoints.size(); ++j){
                    Point i = intersectionPoints.get(j);
                    if (Math.abs(_p.checkPointLocation(i)) < 0.001 && a.isInside(i)){
                        trim = true;
                        if (in == null || Point.distance(in, a.end2) - Point.distance(i, a.end2) > 0.0) {
                            in = i;
                        }
                    }
                }
                if (trim){
                    a.vrts.clear();
                    a.vrts.add(in);
                    a.vrts.add(a.end2);
                    a.setEndPoints(in, a.end2, false);
                    ArcUtil.refineArc(a, Surface.maxEdgeLen, false,0, false);
                    a.baseSubdivision = ArcUtil.getSubdivisionLevel(a);

                    a.valid = true;
                    arc.vrts.clear();
                    arc.vrts.add(arc.end1);
                    arc.vrts.add(in);
                    arc.setEndPoints(arc.end1, in, false);
                    ArcUtil.refineArc(arc, Surface.maxEdgeLen, false,0, false);
                    arc.baseSubdivision = ArcUtil.getSubdivisionLevel(arc);

                    a.prev = arc;
                    arc.next = a;
                    arc.valid = true;
                    newArcs.add(a);
                    newArcs.add(arc);
                    intersectionPoints.remove(in);
                    a.intersecting = false;
                    a = arc;
                    a.intersecting = false;
                } else {
                    a.valid = true;
                    newArcs.add(a);
                }
                a = a.prev;
            } while (a != start);
            for (int i = newArcs.size() - 1; i >= 0; --i){
                b.arcs.add(newArcs.get(i));
            }
            ArcUtil.buildEdges(b, true);
            newBs.add(b);
            sp.boundaries.clear();
            sp.boundaries.addAll(newBs);
        }
        intersectionPoints.clear();
        lastPoI = null;
    }

    public static void processIntersectingConcavePatches(){
        for (SphericalPatch sp : Surface.triangles){
            trimConcavePatch(sp);
        }
    }

    private static void trimConcavePatch(SphericalPatch sp){
        if (!planePoolInitialized){
            for (int i = 0; i < 50; ++i){
                planePool.add(i, new Plane(new Point(0, 0, 0), new Vector(0, 0, 0)));
            }
            planePoolInitialized = true;
        }
        if (vrtsPool.size() == 0){
            for (int i = 0; i < 17; ++i){
                vrtsPool.add(new Point(0, 0, 0));
            }
            _b.arcs.add(_a1);
            _b.arcs.add(_a2);
        }
        try {
            neighbors.clear();
            Surface.probeTree.range(sp.sphere.center.getData(), 2 * SesConfig.probeRadius, neighbors); //causes slf4j warning
            Collections.sort(neighbors, new Comparator<Neighbor<double[], SphericalPatch>>() {
                @Override
                public int compare(Neighbor<double[], SphericalPatch> o1, Neighbor<double[], SphericalPatch> o2) {
                    if (Math.abs(o1.distance - o2.distance) < 0.001) {
                        return 0;
                    }
                    return (o1.distance - o2.distance > 0.0) ? -1 : 1;
                }
            });

            nextPlaneID = 0;
            for (int i = 0; i < sp.boundaries.size(); ++i){
                Boundary b = sp.boundaries.get(i);
                for (int j = 0; j < b.arcs.size(); ++j){
                    Arc a = b.arcs.get(j);
                    if (nextPlaneID >= planePool.size()){
                        planePool.add(new Plane(new Point(0, 0, 0), new Vector(0, 0, 0)));
                    }
                    Plane p = planePool.get(nextPlaneID);
                    p.redefine(a.center, a.normal);
                    nextPlaneID++;
                }
            }
            currIter = 0;
            curr = sp;
            for (int i = 0; i < neighbors.size(); ++i) {
                Neighbor<double[], SphericalPatch> n = neighbors.get(i);
                if (n.value == sp) {
                    continue;
                }
                SphericalPatch sp2 = n.value;
                if (sp.intersectingPatches.contains(sp2.id)) {
                    continue;
                }
                if (nextPlaneID >= planePool.size()){
                    planePool.add(new Plane(new Point(0, 0, 0), new Vector(0, 0, 0)));
                }
                Plane intersectingPlane = planePool.get(nextPlaneID);
                Point center = intersectingPlane.p;
                double radius = computeIntersectionCircle(sp.sphere.center, sp2.sphere.center, center, SesConfig.probeRadius);
                intersectingPlane.v.changeVector(sp.sphere.center, center).makeUnit();
                intersectingPlane.redefine(intersectingPlane.p, intersectingPlane.v);

                if (Point.distance(sp.sphere.center, sp2.sphere.center) < 0.008){
                    continue;
                }
                intersectionPoints.clear();
                exclude.clear();
                findIntersectionPoints(sp, center, radius, intersectionPoints, exclude);

                currInt = intersectionPoints;
                currCirc = intersectingPlane;
                currRad = radius;
                toRemove.clear();
                if (intersectionPoints.size() > 1) {
                    if (!pointsLieOnPatch(sp2, intersectionPoints) && sp.patchNormal.dotProduct(sp2.patchNormal) > 0.8){
                        continue;
                    }
                    boolean identicalIntersector = false;
                    for (int j = 0; j < nextPlaneID; ++j){
                        if (planePool.get(j).isIdenticalWith(intersectingPlane)){
                            identicalIntersector = true;
                            break;
                        }
                    }
                    if (!identicalIntersector){
                        nextPlaneID++;
                    }
                    if (!identicalIntersector){//p!identicalIntersector){//lanes.get(sp.id).stream().noneMatch(plane -> plane.isIdenticalWith(p))) {
                        sp.intersectingPatches.add(sp2.id);

                        generateNewBoundaries(sp, intersectionPoints, intersectingPlane, radius, sp2.id,false);
                        currIter++;
                    }
                } else if (intersectionPoints.size() == 0) {
                    ArcUtil.redefineBoundary(_b, intersectingPlane, radius, vrtsPool, 45);
                    ArcUtil.buildEdges(_b, true);
                    boolean nest = false;
                    removeFromSP.clear();
                    processed.clear();
                    Boundary _newB = null;
                    for (int j = 0; j < sp.boundaries.size(); ++j){
                        Boundary b = sp.boundaries.get(j);
                        if (removeFromSP.contains(b) || processed.contains(b)){
                            continue;
                        }
                        boolean isInside = true;
                        for (int k = 0; k < b.arcs.size(); ++k){
                            Arc a = b.arcs.get(k);
                            rho.redefine(a.center, a.normal);
                            for (int l = 0; l < _b.vrts.size(); ++l){
                                if (!(rho.checkPointLocation(_b.vrts.get(l)) > 0.0)){
                                    isInside = false;
                                    break;
                                }
                            }
                        }
                        if (isInside){
                            for (int k = 0; k < b.nestedBoundaries.size(); ++k){
                                Boundary nb = b.nestedBoundaries.get(k);
                                for (int l = 0; l < nb.arcs.size(); ++l){
                                    Arc a = nb.arcs.get(l);
                                    rho.redefine(a.center, a.normal);
                                    for (int m = 0; m < _b.vrts.size(); ++m){
                                        if (!(rho.checkPointLocation(_b.vrts.get(m)) > 0.0)){
                                            isInside = false;
                                            break;
                                        }
                                    }
                                }
                            }
                            if (isInside) {
                                _newB = new Boundary();
                                Point _center = new Point(_b.arcs.get(0).center);
                                Arc _a1 = ArcUtil.cloneArc(_b.arcs.get(0));
                                _a1.center = _center;
                                Arc _a2 = new Arc(_a1.center, _a1.radius);
                                _a2.setEndPoints(_a1.end2, _a1.end1, false);
                                _a2.setNormal(_a1.normal);
                                _a2.vrts.add(_a2.end1);
                                for (int k = 1; k < _b.arcs.get(1).vrts.size() - 1; ++k){
                                    _a2.vrts.add(new Point(_b.arcs.get(1).vrts.get(k)));
                                }
                                _a2.vrts.add(_a2.end2);
                                _newB.arcs.add(_a1);
                                _newB.arcs.add(_a2);
                                ArcUtil.buildEdges(_newB, true);
                                _newB.nestedBoundaries.add(b);
                                _newB.nestedBoundaries.addAll(b.nestedBoundaries);
                                for (int k = 0; k < _newB.nestedBoundaries.size(); ++k) {
                                    Boundary nb = _newB.nestedBoundaries.get(k);
                                    nb.nestedBoundaries.add(_newB);
                                }
                                toRemove.clear();
                                for (int k = 0; k < _newB.nestedBoundaries.size(); ++k){
                                    Boundary nb = _newB.nestedBoundaries.get(k);
                                    for (int l = 0; l < nb.arcs.size(); ++l){
                                        Arc a = nb.arcs.get(l);
                                        for (int m = 0; m < a.vrts.size(); ++m){
                                            if (!(intersectingPlane.checkPointLocation(a.vrts.get(m)) > 0.0)){
                                                toRemove.add(nb);
                                                break;
                                            }
                                        }
                                    }
                                }
                                _newB.nestedBoundaries.removeAll(toRemove);
                                removeFromSP.addAll(toRemove);
                                for (int k = 0; k < _newB.nestedBoundaries.size(); ++k){
                                    Boundary nb = _newB.nestedBoundaries.get(k);
                                    nb.nestedBoundaries.removeAll(toRemove);
                                }
                                processed.addAll(_newB.nestedBoundaries);
                                nest = true;
                            }
                        }
                    }
                    if (nest){
                        _newB.patch = sp;
                        sp.boundaries.removeAll(removeFromSP);
                        sp.boundaries.add(_newB);
                    }

                }
            }
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private static void trimBoundary(Plane circle, double radius, Point in1, Arc a1, Point in2, Arc a2, List<Point> intersectionPoints, List<Point> usedPoints, List<Boundary> toRemove, List<Boundary> newBS, SphericalPatch sp, int otherPatch){
        Point pStart = in1;
        toRemove2.clear();
        Boundary b = new Boundary();
        Point newArcCenter = new Point(circle.p);
        b.patch = sp;
        Arc start = a1;
        Arc a = a2;

        boolean toBridge = isOptimal(in1, a1, circle);
        boolean forceContinue = false;
        boolean arcExit = false;
        int it = 0;
        do {
            if (toBridge){ //|| forceBridge){
                Arc newA = null; //retrieveBridgeArc(sp.id, otherPatch, in2, in1);
                if (newA == null) {
                    newA = new Arc(newArcCenter, radius);
                    newA.setEndPoints(in1, in2, false);
                    newA.setNormal(circle.v);
                    newA.vrts.add(in1);
                    newA.vrts.add(in2);
                    ArcUtil.refineArc(newA, Surface.maxEdgeLen, false, 0, false);
                } else {
                    Arc temp = newA;
                    newA = ArcUtil.cloneArc(temp);
                    ArcUtil.reverseArc(newA, true);
                    newA.opposite = temp;
                    temp.opposite = newA;
                }
                usedPoints.add(in1);
                usedPoints.add(in2);
                b.arcs.add(newA);
                toBridge = false;
                in1 = in2;
                a1 = a2;
                in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in1, false, circle.p, circle.v, true);
                a2 = in2.arc;//ArcUtil.findContainingArc(in2, circle, sp, null);
                while (a2.bOwner != a1.bOwner && !a2.bOwner.nestedBoundaries.contains(a1.bOwner)){
                    in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in2, false, circle.p, circle.v, true);
                    a2 = ArcUtil.findContainingArc(in2, circle, sp, null);
                }
                a = a1;
                forceContinue = (in1 != pStart);
                arcExit = true;
                if (!toRemove2.contains(a1.bOwner)) {
                    toRemove2.add(a1.bOwner);
                    a1.bOwner.mergeSplit.add(b);
                    b.mergeSplit.add(a1.bOwner);
                }
                if (!toRemove2.contains(a2.bOwner)) {
                    toRemove2.add(a2.bOwner);
                    a2.bOwner.mergeSplit.add(b);
                    b.mergeSplit.add(a2.bOwner);
                }
            } else if (a1 == a2 && ArcUtil.getOrder(a1, in1, in2) < 0) {
                if (ArcUtil.getOrder(a1, in1, in2) < 0) {
                    Arc newA = new Arc(a1.center, a1.radius);
                    newA.setEndPoints(in1, in2, false);
                    newA.setNormal(a1.normal);
                    newA.vrts.add(in1);
                    newA.vrts.add(in2);
                    usedPoints.add(in2);
                    if (a1.opposite != null){
                        newA.opposite = a1.opposite;
                        newA.opposite.opposite = newA;
                        if (a1.torus != null){
                            Arc finalA = a1;
                            a1.torus.concavePatchArcs.removeIf(new Predicate<Arc>() {
                                @Override
                                public boolean test(Arc arc) {
                                    return arc.id == finalA.id;
                                }
                            });
                            newA.torus = a1.torus;
                            newA.torus.concavePatchArcs.add(newA);
                        }
                    }
                    if (a1.cuspTriangle != null){
                        System.out.println("found arc of cusp 1");
                    }

                    ArcUtil.refineArc(newA, Surface.maxEdgeLen, false, 0, false);
                    toBridge = true;
                    in1 = in2;
                    in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in1, false, circle.p, circle.v, true);
                    a2 = in2.arc;//ArcUtil.findContainingArc(in2, circle, sp, null);
                    while (a2.bOwner != a1.bOwner && !a2.bOwner.nestedBoundaries.contains(a1.bOwner)) {
                        in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in2, false, circle.p, circle.v, true);
                        a2 = ArcUtil.findContainingArc(in2, circle, sp, null);
                    }
                    a = a2;
                    b.arcs.add(newA);
                    forceContinue = false;
                    //it is necessary to set a...cause now a==start = true -> premature end
                }
            } else if (a == a1 && arcExit){
                if (Point.distance(in1, a1.end2) > 0.0015) {
                    Arc newA = new Arc(a1.center, a1.radius);
                    newA.setEndPoints(in1, a1.end2, false);
                    newA.setNormal(a.normal);
                    newA.vrts.add(in1);
                    newA.vrts.add(a1.end2);
                    ArcUtil.refineArc(newA, Surface.maxEdgeLen, false, 0, false);
                    if (a.opposite != null){
                        newA.opposite = a.opposite;
                        newA.opposite.opposite = newA;
                        if (a.torus != null){
                            newA.torus = a.torus;
                            Arc finalA = a;
                            a.torus.concavePatchArcs.removeIf(new Predicate<Arc>() {
                                @Override
                                public boolean test(Arc arc) {
                                    return arc.id == finalA.id;
                                }
                            });
                            newA.torus.concavePatchArcs.add(newA);
                        } else if (a.cuspTriangle != null){
                            newA.cuspTriangle = a.cuspTriangle;
                            if (a == a.cuspTriangle.left){
                                newA.cuspTriangle.left = newA;
                                newA.cuspTriangle.right.opposite = newA;
                                newA.opposite = newA.cuspTriangle.right;
                            } else {
                                newA.cuspTriangle.right = newA;
                                newA.cuspTriangle.left.opposite = newA;
                                newA.opposite = newA.cuspTriangle.left;
                            }
                        }
                    }
                    b.arcs.add(newA);
                }
                a = a.next;
                forceContinue = false;
                arcExit = false;
            } else if (a == a2 && !arcExit){
                if (Point.distance(a.end1, in2) > 0.0015) {
                    Arc newA = new Arc(a.center, a.radius);
                    newA.setEndPoints(a.end1, in2, false);
                    newA.setNormal(a.normal);
                    newA.vrts.add(a.end1);
                    newA.vrts.add(in2);
                    ArcUtil.refineArc(newA, Surface.maxEdgeLen, false, 0, false);
                    if (a.opposite != null){
                        newA.opposite = a.opposite;
                        newA.opposite.opposite = newA;
                        if (a.torus != null){//this should not happen
                            newA.torus = a.torus;
                            Arc finalA = a;
                            a.torus.concavePatchArcs.removeIf(new Predicate<Arc>() {
                                @Override
                                public boolean test(Arc arc) {
                                    return arc.id == finalA.id;
                                }

                            });
                            newA.torus.concavePatchArcs.add(newA);
                        } else if (a.cuspTriangle != null){
                            newA.cuspTriangle = a.cuspTriangle;
                            if (a == a.cuspTriangle.left){
                                newA.cuspTriangle.left = newA;
                                newA.cuspTriangle.right.opposite = newA;
                                newA.opposite = newA.cuspTriangle.right;
                            } else {
                                newA.cuspTriangle.right = newA;
                                newA.cuspTriangle.left.opposite = newA;
                                newA.opposite = newA.cuspTriangle.left;
                            }
                        }
                    } else {

                    }
                    b.arcs.add(newA);
                }
                a = a.next;
                forceContinue = false;
                toBridge = true;
                in1 = in2;
                a1 = a2;
                in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in1, false, circle.p, circle.v, true);
                a2 = in2.arc;//ArcUtil.findContainingArc(in2, circle, sp, null);
                while (a2.bOwner != a1.bOwner && !a2.bOwner.nestedBoundaries.contains(a1.bOwner)){
                    in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in2, false, circle.p, circle.v, true);
                    a2 = ArcUtil.findContainingArc(in2, circle, sp, null);
                }
            }else {
                Arc newA = ArcUtil.dbgCloneArc(a);
                if (a.opposite != null){
                    newA.opposite = a.opposite;
                    newA.opposite.opposite = newA;
                    if (a.torus != null){
                        newA.torus = a.torus;
                        int l = 0;
                        while (l < a.torus.concavePatchArcs.size()){
                            if (a.torus.concavePatchArcs.get(l).id == a.id){
                                a.torus.concavePatchArcs.remove(l);
                                break;
                            }
                            l++;
                        }
                        newA.torus.concavePatchArcs.add(newA);
                    } else if (a.cuspTriangle != null){
                        newA.cuspTriangle = a.cuspTriangle;
                        if (a == a.cuspTriangle.left){
                            newA.cuspTriangle.left = newA;
                            newA.cuspTriangle.right.opposite = newA;
                            newA.opposite = newA.cuspTriangle.right;
                        } else {
                            newA.cuspTriangle.right = newA;
                            newA.cuspTriangle.left.opposite = newA;
                            newA.opposite = newA.cuspTriangle.left;
                        }
                    }
                }
                b.arcs.add(newA);
                a = a.next;
            }
        } while (forceContinue || in1 != pStart);
        ArcUtil.buildEdges(b, true);
        newBS.add(b);
        for (int i = 0; i < toRemove2.size(); ++i){
            Boundary b_ = toRemove2.get(i);
            if (!toRemove.contains(b_)){
                toRemove.add(b_);
            }
        }
    }

    private static void generateNewBoundaries(SphericalPatch sp, List<Point> intersectionPs, Plane circle, double radius, int otherPatch, boolean force){
        try {
            if (intersectionPs.size() % 2 == 1){
                return;
            }
            usedPoints.clear();
            toRemove.clear();
            newBS.clear();
            int count = intersectionPoints.size();
            if (intersectionPoints.size() > 1) {
                int i = 0;
                while (i < count) {
                    Point in1 = findOptimalPoint(intersectionPoints, usedPoints, sp, circle);
                    if (in1 == null){
                        break;
                    }
                    Point in2 = ArcUtil.findClosestPointOnCircle(intersectionPoints, in1, false, circle.p, circle.v, true);
                    Point pStart = in1;
                    Arc a1 = in1.arc;//ArcUtil.findContainingArc(in1, circle, sp, null);
                    Arc a2 = in2.arc;//ArcUtil.findContainingArc(in2, circle, sp, null);
                    if (a1 == null || a2 == null){
                        int b = 433;
                    }
                    final Point p1 = in1;
                    final Point p2 = in2;

                    if (a1.bOwner != a2.bOwner && !a1.bOwner.nestedBoundaries.contains(a2.bOwner)){
                        return;
                    }

                    if (Point.distance(a1.end1, p1) < 0.002){
                        if (circle.checkPointLocation(a1.prev.end1) > 0.0){
                            a1 = a1.prev;
                            in1.arc = a1;
                        }
                    }

                    if (Point.distance(a2.end2, p2) < 0.002){
                        if (circle.checkPointLocation(a2.next.end2) > 0.0){
                            a2 = a2.next;
                            in2.arc = a2;
                        }
                    }

                    List<Point> ps = (intersectionPoints.size() > 2) ? getUsablePoints(intersectionPoints, in1, in2, sp, circle) : intersectionPoints;

                    //boolean inside = a1.vrts.stream().allMatch(p -> (Point.distance(p, p1) < 0.001) || circle.checkPointLocation(p) > 0.0 || circle.distanceFromPlane(p) < 0.001) &&
                    //        a1.next.vrts.stream().allMatch(p -> Point.distance(p, p1) < 0.001 || circle.checkPointLocation(p) > 0.0 || circle.distanceFromPlane(p) < 0.001) &&
                    //        a1.prev.vrts.stream().allMatch(p -> Point.distance(p, p1) < 0.001 || circle.checkPointLocation(p) > 0.0 || circle.distanceFromPlane(p) < 0.001) &&
                    //        a2.vrts.stream().allMatch(p -> Point.distance(p, p2) < 0.001 || circle.checkPointLocation(p) > 0.0 || circle.distanceFromPlane(p) < 0.001) &&
                    //        a2.next.vrts.stream().allMatch(p -> Point.distance(p, p2) < 0.001 || circle.checkPointLocation(p) > 0.0 || circle.distanceFromPlane(p) < 0.001) &&
                    //        a2.prev.vrts.stream().allMatch(p -> Point.distance(p, p2) < 0.001 || circle.checkPointLocation(p) > 0.0 || circle.distanceFromPlane(p) < 0.001);
                    boolean inside = false;
                    if (arcInsideHyperspace(a1, circle, p1) && arcInsideHyperspace(a1.next, circle, p1) && arcInsideHyperspace(a1.prev, circle, p1) &&
                            arcInsideHyperspace(a2, circle, p2) && arcInsideHyperspace(a2.next, circle, p2) && arcInsideHyperspace(a2.prev, circle, p2)){
                        inside = true;
                    }
                    if (!force && inside){
                        return;
                    }

                    if (Point.distance(in1, a1.end2) < 0.001 && Point.distance(in2, a2.end2) < 0.001){
                        boolean allInside = true;
                        for (int j = 0; j < a1.bOwner.arcs.size(); ++j){
                            Arc a = a1.bOwner.arcs.get(j);
                            //allInside = allInside && a.vrts.stream().allMatch(p -> (Point.distance(p, p1) < 0.001 || Point.distance(p, p2) < 0.001 || circle.checkPointLocation(p) > 0.0));
                            for (int l = 0; l < a.vrts.size(); ++l){
                                Point point = a.vrts.get(l);
                                if (!(Point.distance(point, p1) < 0.001 || Point.distance(point, p2) < 0.001 || circle.checkPointLocation(point) > 0.0)){
                                    allInside = false;
                                    break;
                                }
                            }
                        }
                        if (allInside){
                            return;
                        }
                    }

                    if (a1.next == a2.prev && Point.distance(in1, a1.end2) < 0.001 && Point.distance(in2, a2.end1) < 0.001){
                        boolean allInside = true;
                        for (int j = 0; j < a1.next.vrts.size(); ++j){
                            Point point = a1.next.vrts.get(j);
                            if (!(Point.distance(point, p1) < 0.001 || Point.distance(point, p2) < 0.001 || circle.checkPointLocation(point) > 0.0)){
                                allInside = false;
                                break;
                            }
                        }
                        if (allInside){
                            return;
                        }
                    }
                    trimBoundary(circle, radius, in1, a1, in2, a2, ps, usedPoints, toRemove, newBS, sp, otherPatch);
                    i += ps.size();
                    intersectionPoints.removeAll(ps);
                }
                updatePatchBoundaries(sp, toRemove, newBS);
                for (int j = 0; j < newBS.size(); ++j){
                    Boundary b = newBS.get(j);
                    if (b.nestedBoundaries.size() > 0){
                        continue;
                    }
                    boolean valid = false;
                    for (int k = 0; k < b.arcs.size(); ++k){
                        Arc a = b.arcs.get(k);
                        if (a.cuspTriangle != null || a.torus != null){
                            valid = true;
                            break;
                        }
                    }
                    if (!valid && SesConfig.verbose){
                        System.out.println("possible invalid boundary for: " + sp.id);
                    }
                }
            }
        } catch (Exception e){
            e.printStackTrace();
            System.err.println("for: " + sp.id);
        }
    }

    private static void findIntersectionPoints(SphericalPatch sp, Point circle, double radius, List<Point> intersectionPoints, List<Arc> exclude){
        circleN.changeVector(sp.sphere.center, circle).makeUnit();
        p1.redefine(circle, circleN);
        for (int i = 0; i < sp.boundaries.size(); ++i){
            Boundary b = sp.boundaries.get(i);
            for (int j = 0; j < b.arcs.size(); ++j){
                Arc a = b.arcs.get(j);
                if (a.torus != null){
                    continue;
                } else if (exclude.size() > 0){
                    boolean _continue = false;
                    for (int k = 0; k < exclude.size(); ++k){
                        Arc a_ = exclude.get(k);
                        if (a_ == a || a_.next == a || a_.prev == a){
                           _continue = true;
                           break;
                        }
                    }
                    if (_continue){
                        continue;
                    }
                }
                p2.redefine(a.center, a.normal);
                if (Plane.getIntersectionLine(p1, p2, vInt, pInt)){
                    hypo.changeVector(circle, pInt);
                    double odvesna = vInt.dotProduct(hypo);
                    double dist = Math.sqrt(hypo.dotProduct(hypo) - odvesna * odvesna);
                    if (dist - radius < 0.0){
                        midOfChord.assignTranslation(pInt, vInt.multiply(odvesna));
                        double odv2 = Math.sqrt(radius * radius - dist * dist);
                        in1.assignTranslation(midOfChord, vInt.makeUnit().multiply(odv2));
                        in2.assignTranslation(midOfChord, vInt.multiply(-1.0));
                        if (a.isInside(in1)){// && ArcUtil.findContainingArc(in1, p1, sp, null) != null){
                            boolean foundSimilar = false;
                            for (int k = 0; k < intersectionPoints.size(); ++k){
                                if (Point.distance(in1, intersectionPoints.get(k)) < 0.01){
                                   foundSimilar = true;
                                   break;
                                }
                            }
                            if (!foundSimilar){//intersectionPoints.stream().noneMatch(v -> Point.distance(in1, v) < 0.01)){
                                Point _p = new Point(in1);
                                _p.arc = a;
                                if (Point.distance(a.end1, _p) < 0.001 && nextSign(_p, a, p1) < 0.0){
                                    _p.arc = a.prev;
                                } else if (Point.distance(a.end2, _p) < 0.001 && nextSign(_p, a, p1) > 0.0){
                                    _p.arc = a.next;
                                }
                                intersectionPoints.add(_p);
                                if (a.cuspTriangle != null){
                                    for (int k = 0; k < a.vrts.size(); ++k){
                                        Point point = a.vrts.get(k);
                                        if (!(p1.checkPointLocation(point) > 0.0 || p1.distanceFromPlane(point) < 0.002)){
                                            intersectionPoints.remove(_p);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if (a.isInside(in2)){// && ArcUtil.findContainingArc(in2, p1, sp, null) != null) {
                            boolean foundSimilar = false;
                            for (int k = 0; k < intersectionPoints.size(); ++k){
                                if (Point.distance(in2, intersectionPoints.get(k)) < 0.01){
                                    foundSimilar = true;
                                    break;
                                }
                            }
                            if (!foundSimilar){//*/intersectionPoints.stream().noneMatch(v -> Point.distance(in2, v) < 0.01)){
                                Point _p = new Point(in2);
                                _p.arc = a;
                                if (Point.distance(a.end1, _p) < 0.001 && nextSign(_p, a, p1) < 0.0){
                                    _p.arc = a.prev;
                                } else if (Point.distance(a.end2, _p) < 0.001 && nextSign(_p, a, p1) > 0.0){
                                    _p.arc = a.next;
                                }
                                intersectionPoints.add(_p);
                                if (a.cuspTriangle != null){
                                    for (int k = 0; k < a.vrts.size(); ++k){
                                        Point point = a.vrts.get(k);
                                        if (!(p1.checkPointLocation(point) > 0.0 || p1.distanceFromPlane(point) < 0.002)){
                                            intersectionPoints.remove(_p);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private static Point findOptimalPoint(List<Point> points, List<Point> usedPoints, SphericalPatch sp, Plane plane){
        try {
            for (int i = 0; i < points.size(); ++i) {
                Point p = points.get(i);
                if (usedPoints.contains(p)){
                    continue;
                }
                Arc a = p.arc;//ArcUtil.findContainingArc(p, plane, sp, null);
                if (a == null){
                    int r = 42;
                }
                if (isOptimal(p, a, plane)){
                    return p;
                }
            }
        } catch (Exception e){
            e.printStackTrace();
        }
        return null;
    }

    private static List<Point> getUsablePoints(List<Point> intersectionPoints, Point origin, Point start, SphericalPatch sp, Plane circle){
        usedPoints2.clear();
        usablePoints.clear();
        usablePoints.addAll(intersectionPoints);
        Arc originArc = origin.arc;//ArcUtil.findContainingArc(origin, circle, sp, null);
        Arc startArc = start.arc;//ArcUtil.findContainingArc(start, circle, sp, null);
        Arc endArc = null;
        Arc a = null;
        Point end = null;
        boolean init = true;
        boolean bridge = (circle.checkPointLocation(startArc.end1) > 0.0);
        do {
            if (init){
                usedPoints2.add(start);
                usablePoints.remove(start);
                startArc = start.arc;//ArcUtil.findContainingArc(start, circle, sp, null);
                end = ArcUtil.findClosestPointOnCircle(intersectionPoints, start, false, circle.p, circle.v, true);
                endArc = end.arc;//ArcUtil.findContainingArc(end, circle, sp, null);
                a = (startArc == endArc) ? ((ArcUtil.getOrder(startArc, start, end) < 0) ? startArc : startArc.next) : startArc.next;
                if (startArc.bOwner != endArc.bOwner && startArc.bOwner.nestedBoundaries.contains(endArc.bOwner)){
                    a = endArc;
                }
                if (end == origin){
                    a = endArc;
                }
                if (isOptimal(start, startArc, circle)){
                    a = endArc;
                }
                init = false;
            } else if (a == endArc){
                start = end;
                init = true;
            } else {
                Point newStart = null;
                for (int i = 0; i < usablePoints.size(); ++i){
                    Point p = usablePoints.get(i);
                    if (a.isInside(p)){
                        Point otherP = null;
                        for (int l = 0; l < usablePoints.size(); ++l){
                            Point v = usablePoints.get(l);
                            if (v != p && a.isInside(v)){
                                otherP = v;
                                break;
                            }
                        }
                        newStart = (otherP != null) ? ((ArcUtil.getOrder(a, p, otherP) < 0) ? p : otherP) : p;
                        break;
                    }
                }
                if (newStart != null){
                    start = newStart;
                    init = true;
                }
                a = a.next;
            }
        } while (init || start != origin);
        return usedPoints2;
    }

    private static void updatePatchBoundaries(SphericalPatch sp, List<Boundary> toRemove, List<Boundary> toAdd){
        for (int i = 0; i < toRemove.size(); ++i){
            Boundary b = toRemove.get(i);
            for (int j = 0; j < b.nestedBoundaries.size(); ++j){
                Boundary nb = b.nestedBoundaries.get(j);
                nb.nestedBoundaries.remove(b);
            }
        }
        for (int i = 0; i < toRemove.size(); ++i){
            Boundary b = toRemove.get(i);
            if (b.mergeSplit.size() == 1){
                Boundary b2 = b.mergeSplit.get(0);
                if (b2.mergeSplit.size() == 1){
                    for (int j = 0; j < b.nestedBoundaries.size(); ++j){
                        Boundary nb = b.nestedBoundaries.get(j);
                        if (areNested(b2, nb)){
                            b2.nestedBoundaries.add(nb);
                            nb.nestedBoundaries.add(b2);
                        } else {
                            sp.boundaries.remove(nb);
                            for (int k = 0; k < nb.nestedBoundaries.size(); ++k){
                                Boundary nb_ = nb.nestedBoundaries.get(k);
                                nb_.nestedBoundaries.remove(nb);
                            }
                        }
                    }
                } else if (b2.mergeSplit.size() > 1){
                    for (int j = 0; j < b.nestedBoundaries.size(); ++j){
                        Boundary nb = b.nestedBoundaries.get(j);
                        if (areNested(b2, nb)){
                            b2.nestedBoundaries.add(nb);
                            nb.nestedBoundaries.add(b2);
                        } else {
                            sp.boundaries.remove(nb);
                            for (int k = 0; k < nb.nestedBoundaries.size(); ++k){
                                Boundary nb_ = nb.nestedBoundaries.get(k);
                                nb_.nestedBoundaries.remove(nb);
                            }
                        }
                    }

                }
                for (int j = 0; j < b2.mergeSplit.size(); ++j){
                    b2.mergeSplit.get(j).mergeSplit.clear();
                }
            } else if (b.mergeSplit.size() > 1){
                for (int j = 0; j < b.mergeSplit.size(); ++j){
                    Boundary newB = b.mergeSplit.get(j);
                    for (int k = 0; k < b.nestedBoundaries.size(); ++k){
                        Boundary nb = b.nestedBoundaries.get(k);
                        if (areNested(newB, nb)){
                            newB.nestedBoundaries.add(nb);
                            nb.nestedBoundaries.add(newB);
                        } else {
                            sp.boundaries.remove(nb);
                            for (int l = 0; l < nb.nestedBoundaries.size(); ++l){
                                Boundary nb_ = nb.nestedBoundaries.get(l);
                                nb_.nestedBoundaries.remove(nb);
                            }
                        }
                    }
                    newB.mergeSplit.clear();
                }
                b.mergeSplit.clear();
            }
        }
        for (int i = 0; i < toAdd.size(); ++i){
            toAdd.get(i).mergeSplit.clear();
        }
        sp.boundaries.removeAll(toRemove);
        sp.boundaries.addAll(toAdd);
    }

    private static boolean areNested(Boundary b1, Boundary b2){
        boolean _nest = true;
        for (int i = 0; i < b1.arcs.size(); ++i){
            Arc a = b1.arcs.get(i);
            rho.redefine(a.center, a.normal);
            for (int j = 0; j < b2.vrts.size(); ++j){
                if (!(rho.checkPointLocation(b2.vrts.get(j)) > 0.0)){
                    _nest = false;
                    break;
                }
            }
            if (!_nest){
                break;
            }
        }
        if (_nest){
            for (int i = 0; i < b2.arcs.size(); ++i){
                Arc a = b2.arcs.get(i);
                rho.redefine(a.center, a.normal);
                for (int j = 0; j < b1.vrts.size(); ++j){
                    if (!(rho.checkPointLocation(b1.vrts.get(j)) > 0.0)){
                        _nest = false;
                        break;
                    }
                }
                if (!_nest){
                    break;
                }
            }
            if (!_nest){
                return false;
            }
        } else {
            return false;
        }
        return true;
    }

    private static boolean pointsLieOnPatch(SphericalPatch sp, List<Point> points){
        for (int i = 0; i < sp.boundaries.size(); ++i){
            Boundary b = sp.boundaries.get(i);
            for (int j = 0; j < b.arcs.size(); ++j){
                Arc a = b.arcs.get(j);
                for (int k = 0; k < points.size(); ++k){
                    Point p = points.get(k);
                    if (a.isInside(p)){
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static boolean isOptimal(Point p, Arc a, Plane circle){
        int dir = 1;
        if (Point.distance(p, a.end2) < 0.001){
            dir = -1;
        }
        Point _p = genP(a, p, dir);
        _genVector.changeVector(_nextPoint, circle.p).multiply(10.f);
        _nextPoint.assignTranslation(circle.p, _genVector);
        return dir * circle.checkPointLocation(_nextPoint) < 0.0;
    }

    public static double nextSign(Point p, Arc a, Plane circle){
        Point _p = genP(a, p, 1);
        _genVector.changeVector(_p, circle.p).multiply(10.f);
        _nextPoint.assignTranslation(circle.p, _genVector);
        return circle.checkPointLocation(_nextPoint);
    }

    public static Point genP(Arc a, Point p, int dir){
        _genVector.changeVector(p, a.center).makeUnit();
        _genQuaternion.setFromAngleNormalAxis((float)(dir * Math.toRadians(2)), a.normal.getFloatData());
        _genQuaternion.rotateVector(_floatVector, 0, _genVector.getFloatData(), 0);
        _genVector.changeVector(_floatVector[0], _floatVector[1], _floatVector[2]);
        _genVector.makeUnit().multiply(a.radius);
        return _nextPoint.assignTranslation(a.center, _genVector);
    }

    private static boolean arcInsideHyperspace(Arc a, Plane circle, Point p){
        boolean inside = true;
        for (int i = 0; i < a.vrts.size(); ++i){
            Point point = a.vrts.get(i);
            if (!(Point.distance(p, point) < 0.001 || circle.checkPointLocation(point) > 0.0 || circle.distanceFromPlane(point) < 0.001)){
               inside = false;
               break;
            }
        }
        return inside;
    }

    private static void refineArcForNvertices(Arc a, int n){
        if (a.vrts.size() >= n){
            return;
        }
        int currSubd = ArcUtil.getSubdivisionLevel(a);
        int target = (int)(Math.log10(n) / Math.log10(2));
        ArcUtil.refineArc(a, SesConfig.edgeLimit, true, (target - currSubd), false);
    }

    private static void equalizeVertexCount(CuspTriangle tr1, CuspTriangle tr2){
        int maxCount = Math.max(Math.max(tr1.left.vrts.size(), tr1.right.vrts.size()), Math.max(tr2.left.vrts.size(), tr2.right.vrts.size()));
        refineArcForNvertices(tr1.left, maxCount);
        refineArcForNvertices(tr1.right, maxCount);
        refineArcForNvertices(tr2.left, maxCount);
        refineArcForNvertices(tr2.right, maxCount);
    }

    public static int getTorusProbeIdx(ToroidalPatch tp, int vertexID){
        if (tp.tr1 != null){
           int firstPatchVertexCount = tp.tr1.base.vrts.size() * tp.tr1.left.vrts.size();
           int arcLen = (vertexID >= firstPatchVertexCount) ? tp.tr2.left.vrts.size() : tp.tr1.left.vrts.size();
           int probeOffset = (vertexID >= firstPatchVertexCount) ? tp.tr1.base.vrts.size() : 0;
           int vertexOffset = (probeOffset > 0) ? firstPatchVertexCount : 0;
           return (vertexID - vertexOffset) / arcLen + probeOffset;
        } else {
           if (tp.concavePatchArcs.size() == 0){
              if (getProbeAxisDistance(tp.probe1, tp.convexPatchArcs.get(0).center, tp.convexPatchArcs.get(1).center) - SesConfig.probeRadius < 0.0){
                  int arcLen = tp.vertices.size() / (2 * tp.convexPatchArcs.get(0).vrts.size());
                  return vertexID / arcLen;
              } else {
                 int arcLen = tp.vertices.size() / tp.convexPatchArcs.get(0).vrts.size();
                 return vertexID / arcLen;
              }
           } else {
               return vertexID / tp.concavePatchArcs.get(0).vrts.size();
           }
        }
    }

    public static double computeIntersectionCircle(Point probe1, Point probe2, Point result, double probeRadius){
        v1.changeVector(probe1, probe2).multiply(0.5f);
        point.assignTranslation(probe2, v1);
        result.x = point.x;
        result.y = point.y;
        result.z = point.z;
        return Math.sqrt(Math.pow(probeRadius, 2) - Math.pow(v1.sqrtMagnitude(), 2));
    }

    /*
    computes cusp points for self-intersecting toroidal patch
     */
    public static Point computeCusp(Point probe, Sphere a1, Sphere a2){
        v1.changeVector(a2.center, a1.center).makeUnit(); //axis
        v2.changeVector(probe, a1.center).makeUnit(); //a1toprobe
        double atomToProbeLength = a1.radius + SesConfig.probeRadius;
        v2.makeUnit();
        double alpha = Math.acos(v1.dotProduct(v2));
        double beta = Math.asin(((atomToProbeLength) * Math.sin(alpha)) / SesConfig.probeRadius);
        double gama = Math.PI - alpha - beta;
        double atomToCusp = (SesConfig.probeRadius * Math.sin(gama)) / Math.sin(alpha);
        v1.multiply(atomToCusp);
        return Point.translatePoint(a1.center, v1);
    }

    public static double getProbeAxisDistance(Point probe, Point a1, Point a2){
        v1.changeVector(a2, a1).makeUnit(); //axis vector
        v2.changeVector(probe, a1); //vector to probe
        v1.multiply(v2.dotProduct(v1));
        return v3.assignAddition(v2, v1.multiply(-1)).sqrtMagnitude();
    }
}

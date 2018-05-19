package cz.fi.muni.xmraz3.mesh;


import com.jogamp.opengl.math.Quaternion;
import cz.fi.muni.xmraz3.SesConfig;
import cz.fi.muni.xmraz3.Surface;
import cz.fi.muni.xmraz3.math.Plane;
import cz.fi.muni.xmraz3.math.Point;
import cz.fi.muni.xmraz3.math.Vector;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class AdvancingFrontMethod {

    private boolean patchComplete = true;
    private boolean atomComplete = true;
    private boolean boundaryComplete = true;
    private boolean loopDetected = false;
    private boolean faceGenerated = false;
    private boolean verbose = false;
    public boolean loop = false;

    private static int timeout = 5000;
    private static int maxNumberOfRestarts = 0;
    private int currentTry = 0;
    private int activeLoop = 0;
    private int numOfLoops = 1;
    private int nextEdgeID = 0;
    private int nextFaceID = 0;
    private int nextNEMID = 0;
    private int numOfTriangles = 0;

    private double edgeLength;
    private double distTolerance;
    private double height;
    private double pointEdgeDistTolerance;
    private double minAlpha;

    private Vector aV1 = new Vector(0, 0, 0);
    private Vector aV2 = new Vector(0, 0, 0);
    private Vector v3 = new Vector(0, 0, 0);
    private Vector n = new Vector(0, 0, 0);
    private Vector ve1 = new Vector(0, 0, 0);
    private Vector ve2 = new Vector(0, 0, 0);
    private Vector v2e1 = new Vector(0, 0, 0);
    private Vector v2e2 = new Vector(0, 0, 0);
    private Vector addv1 = new Vector(0, 0, 0);
    private Vector addv2 = new Vector(0, 0, 0);
    private Vector intersection = new Vector(0, 0, 0);
    private Vector n1 = new Vector(0, 0, 0);
    private Vector n2 = new Vector(0, 0, 0);
    private Vector e1ToE2 = new Vector(0, 0, 0);
    private Vector e2ToE1 = new Vector(0, 0, 0);
    private Vector midVec = new Vector(0, 0, 0);
    private Vector midNormal = new Vector(0, 0, 0);
    private Vector tau1 = new Vector(0, 0, 0);
    private Vector tau2 = new Vector(0, 0, 0);
    private Vector midEF1 = new Vector(0, 0, 0);
    private Vector midEF2 = new Vector(0, 0, 0);
    private Vector tangentInMiddle = new Vector(0, 0, 0);

    private Point _mid1 = new Point(0, 0, 0);
    private Point _mid2 = new Point(0, 0, 0);
    private Point midPoint = new Point(0, 0, 0);
    private Point testPoint = new Point(0, 0, 0);

    private Quaternion q = new Quaternion();
    private float[] nvector = new float[3];
    private float[] _vecToRotate = new float[3];

    private Plane ro1 = new Plane(new Point(0, 0, 0), new Vector(0, 0, 0));
    private Plane ro2 = new Plane(new Point(0, 0, 0), new Vector(0, 0, 0));

    private Edge e = null;
    private Edge e1 = new Edge(0, 0);
    private Edge e2 = new Edge(0, 0);
    private Edge pfp1 = new Edge(0, 0);
    private Edge pfp2 = new Edge(0, 0);
    private Edge ef1E1 = new Edge(0, 0);
    private Edge ef1E2 = new Edge(0, 0);
    private Edge ef2E1 = new Edge(0, 0);
    private Edge ef2E2 = new Edge(0, 0);

    private List<Point> removePoints = new ArrayList<>();
    private List<Edge> facets;
    private List<Edge> pastFacets;
    private List<Point> nodes;
    private List<List<Edge>> nodeEdgeMap;
    private List<Edge> ignore = new ArrayList<>();
    private List<Boundary> processedBoundaries;
    private List<Point> candidates = new ArrayList<>();
    private List<Boundary> toProcess = new ArrayList<>();
    private List<Edge> edgePool = new ArrayList<>(500);
    private List<Face> facePool = new ArrayList<>(300);
    private List<Edge> relevantEdges = new ArrayList<>();
    private List<Edge> boundaryEdges = new ArrayList<>();
    private List<Point> trueCands = new ArrayList<>();
    private LinkedList<Edge> loops;

    private SphericalPatch patch;

    public AdvancingFrontMethod(){
        facets = new ArrayList<>();
        nodes = new ArrayList<>();
        pastFacets = new ArrayList<>();
        nodeEdgeMap = new ArrayList<>(150);
        processedBoundaries = new ArrayList<>();
        loops = new LinkedList<>();
        for (int i = 0; i < 500; ++i){
            edgePool.add(i, new Edge(0, 0));
        }
        for (int i = 0; i < 300; ++i){
            facePool.add(i, new Face(0, 0, 0));
        }
        for (int i = 0; i < 150; ++i){
            nodeEdgeMap.add(new ArrayList<>());
        }
    }

    /*
    One time initialization per patch
     */
    private void initializeAFM(SphericalPatch sp){
        edgeLength = SesConfig.edgeLimit;
        loopDetected = false;
        minAlpha = SesConfig.minAlpha;
        distTolerance = 0.2 * SesConfig.edgeLimit;
        pointEdgeDistTolerance = 0.4 * SesConfig.edgeLimit;
        height = SesConfig.edgeLimit * (Math.sqrt(3) / 2);
        verbose = false;
        facets.clear();
        nodes.clear();
        pastFacets.clear();
        nextEdgeID = 0;
        nextNEMID = 0;
        nextFaceID = 0;
        patch = sp;
        processedBoundaries.clear();
        patchComplete = false;
        boundaryComplete = true;
    }

    /*
    Initalizes data structures. Executes several times depending on the number of independent boundaries(not nested) in patch
     */
    private void initializeAFMDataStructures(){
        //if (atomComplete){
        //    processedBoundaries.clear();
        //    atomComplete = false;
        //}
        loopDetected = false;
        facets.clear();
        pastFacets.clear();
        nodes.clear();
        loops.clear();
        nextEdgeID = 0;
        nextNEMID = 0;
        Boundary b = null;
        for (int i = 0; i < patch.boundaries.size(); ++i){
            Boundary c = patch.boundaries.get(i);
            if (!processedBoundaries.contains(c)){
                b = c;
                break;
            }
        }
        int insideBoundaryCount = 0;
        insideBoundaryCount = b.nestedBoundaries.size();
        for (int i = 0; i < patch.boundaries.size(); ++i){
            Boundary c = patch.boundaries.get(i);
            if (processedBoundaries.contains(c)){
                continue;
            }
            if (c.nestedBoundaries.size() >= insideBoundaryCount){
                insideBoundaryCount = c.nestedBoundaries.size();
                b = c;
            }
        }
        processedBoundaries.add(b);
        processedBoundaries.addAll(b.nestedBoundaries);
        toProcess.clear();
        toProcess.add(b);
        toProcess.addAll(b.nestedBoundaries);
        for (int i = 0; i < toProcess.size(); ++i) {
            Boundary c = toProcess.get(i);
            boundaryEdges.clear();
            for (int j = 0; j < c.vrts.size(); ++j){
                Point p = c.vrts.get(j);
                if (nextNEMID >= nodeEdgeMap.size()){
                    nodeEdgeMap.add(new ArrayList<>());
                }
                p.afmIdx = nextNEMID;
                nodeEdgeMap.get(nextNEMID).clear();
                nextNEMID++;
            }
            for (int j = 0; j < c.vrts.size(); ++j){
                Point p1 = c.vrts.get(j);
                Point p2 = (j < c.vrts.size() - 1) ? c.vrts.get(j + 1) : c.vrts.get(0);
                if (nextEdgeID >= edgePool.size()){
                    edgePool.add(new Edge(0, 0));
                }
                Edge e = edgePool.get(nextEdgeID++);
                e.p1 = p1;
                e.p2 = p2;
                e.v1 = p1.afmIdx;
                e.v2 = p2.afmIdx;
                nodeEdgeMap.get(p1.afmIdx).add(e);
                nodeEdgeMap.get(p2.afmIdx).add(e);
                boundaryEdges.add(e);
            }
            nodes.addAll(c.vrts);
            for (int j = 0; j < boundaryEdges.size(); ++j){
                Edge e1 = boundaryEdges.get(j);
                Edge e2 = (j < boundaryEdges.size() - 1) ? boundaryEdges.get(j + 1) : boundaryEdges.get(0);
                e1.next = e2;
                e2.prev = e1;
            }
            facets.addAll(boundaryEdges);
        }
        e = facets.get(0);
        activeLoop = 0;
        numOfLoops = 1;
        patchComplete = false;
        boundaryComplete = false;
    }

    public void meshSphericalPatch(SphericalPatch sp){
        initializeAFM(sp);
        do {
            this.generateMesh();
        } while (!patchComplete);
    }


    private boolean generateMesh(){
        loop = false;
        if (boundaryComplete){
            this.initializeAFMDataStructures();
        }
        Long time = System.currentTimeMillis();
        int empty = 0;
        while (facets.size() > 0){
            pastFacets.clear();
            faceGenerated = false;
            if (empty >= facets.size()){
                System.out.println("DETECED LOOP, ending mesh generation for patch id: " + patch.id);
                patchComplete = true;
                loopDetected = true;
                break;
            }
            //if (System.currentTimeMillis() - time > timeout){
            //    if (currentTry < maxNumberOfRestarts){
            //        time = System.currentTimeMillis();
            //        currentTry++;
            //    } else {
            //        System.out.println("DETECTED LOOP, ending mesh generation");
            //        atomComplete = true;
            //        loopDetected = true;
            //        break;
            //    }
            //}
            candidates.clear();
            if (e.next.next == e.prev){//only three edges in the current loop -> close it with one triangle
                closeLoop();
                if (loops.size() == 0 && facets.size() > 0){
                    e = facets.get(0);
                } else {
                    while (!facets.contains(e) && loops.size() > 0) {
                        e = loops.poll();
                    }
                }
                activeLoop = e.loopID;
                if (facets.size() == 0){
                    break;
                }
                continue;
            }
            activeLoop = e.loopID;

            n1 = n1.changeVector(e.p1, patch.sphere.center).makeUnit();
            n2 = n2.changeVector(e.p2, patch.sphere.center).makeUnit();
            e1ToE2 = e1ToE2.changeVector(e.p2, e.p1);
            e2ToE1 = e2ToE1.changeVector(e.p1, e.p2);
            midVec = midVec.changeVector(e.p2, e.p1).multiply(0.5);
            midPoint.assignTranslation(e.p1, midVec);
            midNormal = midNormal.changeVector(midPoint, patch.sphere.center).makeUnit();
            tangentInMiddle.assignNormalVectorOf(midNormal, midVec).makeUnit();
            double realHeight = Math.sqrt(Math.pow(edgeLength, 2) - Math.pow(e1ToE2.sqrtMagnitude() * 0.5, 2));
            double realMinAlpha = Math.asin(realHeight / edgeLength) + Math.toRadians(1);
            realMinAlpha = (realMinAlpha > Math.toRadians(120)) ? Math.toRadians(120) : realMinAlpha;

            Edge eR = e.next;
            Edge eL = e.prev;
            //compute the angle between the edge e and e.next, and the angle between the edge e and e.prev
            double alpha1 = computeAngle(aV1.changeVector(e.p1, e.p2).makeUnit(), aV2.changeVector(e.next.p2, e.next.p1).makeUnit(), n2);
            double alpha2 = computeAngle(aV1.changeVector(e.prev.p1, e.prev.p2).makeUnit(), aV2.changeVector(e.p2, e.p1).makeUnit(), n1);

            if (alpha1 < realMinAlpha) {//realMinAlpha){//this.minAlpha){
                if (nodeEdgeMap.get(eR.p2.afmIdx).size() > 0) {
                    e1.p1 = e.p1;
                    e1.p2 = eR.p2;
                    e2.p1 = e.p2;
                    e2.p2 = eR.p2;
                    ignore.clear();
                    ignore.add(e);
                    boolean intersects = false;
                    for (int i = 0; i < facets.size(); ++i){
                        Edge ek = facets.get(i);
                        if (ek == e || ek == e.prev || ek == e.next || ek.loopID != activeLoop){
                            continue;
                        }
                        if (checkForIntersectingEdges(e1, ek, patch.sphere.radius, patch.sphere.center)){
                            if (pointEdgeDistance(ek.p1, e1) - pointEdgeDistance(ek.p2, e1) > 0.0){
                                candidates.add(ek.p2);
                            } else {
                                candidates.add(ek.p1);
                            }
                            intersects = true;
                            break;
                        }
                    }
                    if (!intersects){
                        candidates.add(eR.p2);
                    }
                }
            }
            if (alpha2 < realMinAlpha){
                if (nodeEdgeMap.get(eL.p1.afmIdx).size() > 0) {
                    e1.p1 = e.p1;
                    e1.p2 = eL.p1;
                    e2.p1 = e.p2;
                    e2.p2 = eL.p1;
                    ignore.clear();
                    ignore.add(e);
                    boolean intersects = false;
                    for (int i = 0; i < facets.size(); ++i){
                        Edge ek = facets.get(i);
                        if (ek == e || ek == e.prev || ek == e.next || ek.loopID != activeLoop){
                            continue;
                        }
                        if (checkForIntersectingEdges(e2, ek, patch.sphere.radius, patch.sphere.center)){
                            if (pointEdgeDistance(ek.p1, e2) - pointEdgeDistance(ek.p2, e2) > 0.0){
                                candidates.add(ek.p2);
                            } else {
                                candidates.add(ek.p1);
                            }
                            intersects = true;
                            break;
                        }
                    }
                    if (!intersects){
                        candidates.add(eL.p1);
                    }
                }
            }

            trueCands.clear();
            if (candidates.isEmpty() || true) {
                criterionPointEdgeDistance();
            }
            while (candidates.contains(e.p1)) {
                candidates.remove(e.p1);
            }
            while (candidates.contains(e.p2)){
                candidates.remove(e.p2);
            }
            for (int i = 0; i < candidates.size(); ++i){
                Point p = candidates.get(i);
                if (computeAngle(aV1.changeVector(p, e.p1).makeUnit(), aV2.changeVector(e.p2, e.p1).makeUnit(), n1) > Math.toRadians(SesConfig.minAlpha)
                        || computeAngle(aV1.changeVector(e.p1, e.p2).makeUnit(), aV2.changeVector(p, e.p2).makeUnit(), n2) > Math.toRadians(SesConfig.minAlpha)){// || (Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, p))) < 0.0)) {
                    continue;
                }
                trueCands.add(p);
            }
            candidates.clear();
            candidates.addAll(trueCands);
            //so far no suitable points to form a new triangle with the edge e, let's create a test point
            if (candidates.isEmpty()){
                double height2; //assign such height so that the newly formed edges will have +- length of Surface.maxEdgeLen
                if (edgeLength - 0.05f - e1ToE2.sqrtMagnitude() * 0.5f < 0.f){
                    height2 = (e1ToE2.sqrtMagnitude() * 0.5f) * Math.tan(Math.toRadians(30));
                } else {
                    height2 = Math.sqrt(Math.pow(edgeLength - 0.05f, 2) - Math.pow(e1ToE2.sqrtMagnitude() * 0.5f, 2));
                }
                Point pTest = generateNewTestPoint(midNormal, e1ToE2, patch.sphere.radius, height2, midPoint, true);
                //test whether there are any edges close enough to pTest so the triangle might be formed with them instead of pTest
                pointEdgeDistanceCriterion(pTest);
                //test whether there are any points close enough to pTest so they can be used to form the new triangle
                pointPointDistanceCriterion(pTest);
                while (candidates.contains(e.p1)) {
                    candidates.remove(e.p1);
                }
                while (candidates.contains(e.p2)){
                    candidates.remove(e.p2);
                }
                for (int i = 0; i < candidates.size(); ++i){
                    Point c = candidates.get(i);
                    if (Math.abs(aV1.changeVector(c, e.p1).makeUnit().dotProduct(aV2.changeVector(c, e.p2).makeUnit()) - 1.0) < 0.01){
                        removePoints.add(c);
                    }
                }
                candidates.removeAll(removePoints);
                removePoints.clear();
                if (candidates.isEmpty()) {
                    e1.p1 = e.p1;
                    e1.p2 = pTest;
                    e2.p1 = e.p2;
                    e2.p2 = pTest;
                    ignore.clear();
                    ignore.add(e);
                    if (checkForIntersectingEdges(e1, e2, facets, ignore)){// || checkForIntersectingEdges(e1, e2, pastFacets, ignore)) {
                        if (!(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, e.next.p2))) < 0.01) && hasCorrectOrientation(e.p1, e.p2, e.next.p2)){
                            boolean intersects = false;
                            for (int i = 0; i < facets.size(); ++i){
                                Edge ek = facets.get(i);
                                if (ek == e || ek == e.prev || ek == e.next){// || ek.loopID != activeLoop){
                                    continue;
                                }
                                if (checkForIntersectingEdges(eR, ek, patch.sphere.radius, patch.sphere.center)){
                                    if (pointEdgeDistance(ek.p1, eR) - pointEdgeDistance(ek.p2, eR) > 0.0){
                                        candidates.add(ek.p2);
                                    } else {
                                        candidates.add(ek.p1);
                                    }
                                    intersects = true;
                                    break;
                                }
                            }
                            if (!intersects){
                                candidates.add(e.next.p2);
                                if (verbose){
                                    System.out.println("Added next edge to candidates later");
                                }
                            }
                        }
                        if (!(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, e.prev.p1))) < 0.01) && hasCorrectOrientation(e.p1, e.p2, e.prev.p1)){
                            boolean intersects = false;
                            for (int i = 0; i < facets.size(); ++i){
                                Edge ek = facets.get(i);
                                if (ek == e || ek == e.prev || ek == e.next){// || ek.loopID != activeLoop){
                                    continue;
                                }
                                if (checkForIntersectingEdges(eL, ek, patch.sphere.radius, patch.sphere.center)){
                                    if (pointEdgeDistance(ek.p1, eL) - pointEdgeDistance(ek.p2, eL) > 0.0){
                                        candidates.add(ek.p2);
                                    } else {
                                        candidates.add(ek.p1);
                                    }
                                    intersects = true;
                                    break;
                                }
                            }
                            if (!intersects){
                                candidates.add(e.prev.p1);
                                if (verbose){
                                    System.out.println("Added prev edge to candidates later");
                                }
                            }
                        }
                    } else {
                        generateFaceWithNewPoint();
                        candidates.clear();
                        faceGenerated = true;
                    }
                }
            }
            if (!candidates.isEmpty()){
                Point pt = candidates.get(0);
                double min = pointEdgeDistance(pt, e);
                for (int i = 1; i < candidates.size(); ++i){
                    Point t = candidates.get(i);
                    if (t == e.p1 || t == e.p2)
                        continue;
                    if (pointEdgeDistance(t, e) < min && nodeEdgeMap.get(t.afmIdx).size() > 0){
                        min = pointEdgeDistance(t, e);
                        pt = t;
                    }
                }
                if (pt == e.prev.p1){
                    generateFaceWithPreviousEdge();
                } else if (pt == e.next.p2){
                    generateFaceWithNextEdge();
                } else {
                    generateBridgeFace(pt);
                }
            } else {

            }
            if (!faceGenerated){
                e = e.next;
                empty++;
            } else {
                empty = 0;
            }
        }
        if (facets.size() == 0){
            boundaryComplete = true;
            if (processedBoundaries.size() == patch.boundaries.size()){
                patchComplete = true;
            }
        }
        if (loopDetected){
            patchComplete = true;
            atomComplete = true;
            loop = true;
        }
        return patchComplete;
    }

    public void transferFacesToPatch(){
        patch.faces = new int[3 * nextFaceID];
        int j = 0;
        for (int i = 0; i < nextFaceID; ++i){
            Face f = facePool.get(i);
            patch.faces[j] = f.a;
            patch.faces[j + 1] = f.b;
            patch.faces[j + 2] = f.c;
            j += 3;
        }
        nextFaceID = 0;
    }


    private void criterionPointEdgeDistance(){
        for (int i = 0; i < facets.size(); ++i){
            Edge ef = facets.get(i);
            if (ef == e){ // || ef.loopID != activeLoop){
                continue;
            }
            if (ef == e.prev || ef == e.next){
                continue;
            }
            if (/*ef.p1 != e.p1 && ef.p1 != e.p2*/true){
                double dist = pointEdgeDistance(ef.p1, e);
                if (dist < 2 * edgeLength + pointEdgeDistTolerance){
                    //Vector midToefP1 = Point.subtractPoints(ef.p1, midPoint).makeUnit();
                    aV1.changeVector(ef.p1, midPoint).makeUnit();
                    if (aV1.dotProduct(tangentInMiddle) > 0.0){
                        pfp1.p1 = ef.p1;
                        pfp1.p2 = e.p1;
                        pfp2.p1 = ef.p1;
                        pfp2.p2 = e.p2;
                        ignore.clear();
                        ignore.add(e);
                        ignore.add(ef);
                        if (!checkForIntersectingEdges(pfp1, pfp2, facets, ignore) && !checkForIntersectingEdges(pfp1, pfp2, pastFacets, ignore)){
                            if (nodeEdgeMap.get(ef.p1.afmIdx).size() > 0 && !(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, ef.p1))) < 0.01)) {
                                if (ef.p1 != e.p2) {
                                    candidates.add(ef.p1);
                                    //ef.p1.afmSelect = 2;
                                }
                            }
                        }
                    }
                }
            }
            if (/*ef.p2 != e.p1 && ef.p2 != e.p2*/true){
                double dist = pointEdgeDistance(ef.p2, e);
                if (dist < 2 * edgeLength + pointEdgeDistTolerance){
                    //Vector midToefP1 = Point.subtractPoints(ef.p2, midPoint).makeUnit();
                    aV1.changeVector(ef.p2, midPoint).makeUnit();
                    if (aV2.dotProduct(tangentInMiddle) > 0){
                        pfp1.p1 = ef.p2;
                        pfp1.p2 = e.p1;
                        pfp2.p1 = ef.p2;
                        pfp2.p2 = e.p2;
                        ignore.clear();
                        ignore.add(e);
                        ignore.add(ef);
                        if (!checkForIntersectingEdges(pfp1, pfp2, facets, ignore) && !checkForIntersectingEdges(pfp1, pfp2, pastFacets, ignore)){
                            if (nodeEdgeMap.get(ef.p2.afmIdx).size() > 0 && !(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, ef.p2))) < 0.01)) {
                                if (ef.p2 != e.p1) {
                                    candidates.add(ef.p2);
                                    //ef.p2.afmSelect = 2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void pointEdgeDistanceCriterion(Point p){
        for (int i = 0; i < facets.size(); ++i){
            Edge eF = facets.get(i);
            if (eF == e || eF.loopID != activeLoop){// || eF.loopID != activeLoop){
                continue;
            }
            double pEdgeDist = pointEdgeDistance(p, eF);
            if (Point.distance(eF.p2, eF.p1) + pointEdgeDistTolerance > pEdgeDist){
                midEF1 = midEF1.changeVector(eF.p1, midPoint);
                midEF2 = midEF2.changeVector(eF.p2, midPoint);
                if (midEF1.dotProduct(tangentInMiddle) > 0 && midEF2.dotProduct(tangentInMiddle) > 0){
                    ef1E1.p1 = eF.p1;
                    ef1E1.p2 = e.p1;
                    ef1E2.p1 = eF.p1;
                    ef1E2.p2 = e.p2;
                    ef2E1.p1 = eF.p2;
                    ef2E1.p2 = e.p1;
                    ef2E2.p1 = eF.p2;
                    ef2E2.p2 = e.p2;
                    ignore.clear();
                    ignore.add(e);
                    ignore.add(eF);
                    if (!checkForIntersectingEdges(ef1E1, ef1E2, pastFacets, ignore) && !checkForIntersectingEdges(ef1E1, ef1E2, facets, ignore) && nodeEdgeMap.get(eF.p1.afmIdx).size() > 0 && eF.p1 != e.p2 && !(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, eF.p1))) < 0.01)){
                        candidates.add(eF.p1);
                    }
                    if (!checkForIntersectingEdges(ef2E1, ef2E2, pastFacets, ignore) && !checkForIntersectingEdges(ef2E1, ef2E2, facets, ignore) && nodeEdgeMap.get(eF.p2.afmIdx).size() > 0 && eF.p2 != e.p1  && !(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, eF.p2))) < 0.01)){
                        candidates.add(eF.p2);
                    }
                }
            }
        }
    }

    private void pointPointDistanceCriterion(Point p){
        for (int i = 0; i < nodes.size(); ++i){
            Point p2 = nodes.get(i);
            if (p2 == e.p1 || p2 == e.p2){
                continue;
            }
            double dist = Point.distance(p, p2);
            if (dist < this.distTolerance){
                pfp1.p1 = p2;
                pfp1.p2 = e.p1;
                pfp2.p1 = p2;
                pfp2.p2 = e.p2;
                ignore.clear();
                ignore.add(e);
                if (!checkForIntersectingEdges(pfp1, pfp2, facets, ignore) && !checkForIntersectingEdges(pfp1, pfp2, pastFacets, ignore)){
                    if (nodeEdgeMap.get(p2.afmIdx).size() > 0 && !(Math.abs(midNormal.dotProduct(computeTriangleNormal(e.p1, e.p2, p2))) < 0.01)) {
                        candidates.add(p2);
                    }
                }
            }
        }
    }


    private Point generateNewTestPoint(Vector normal, Vector edgeVector, double sphereRadius, double height, Point origin, boolean t){
        double cosA = (2 * Math.pow(sphereRadius, 2) - Math.pow(height, 2)) / (2 * Math.pow(sphereRadius, 2));
        double angle = Math.PI - Math.acos(cosA);
        angle *= 0.5f;
        angle = Math.PI - angle;
        aV1.changeVector(edgeVector).makeUnit();
        _vecToRotate[0] = (float)aV1.getX();
        _vecToRotate[1] = (float)aV1.getY();
        _vecToRotate[2] = (float)aV1.getZ();
        q.setFromAngleNormalAxis(-1.0f * (float)angle, _vecToRotate);
        nvector[0] = (float)normal.getX();
        nvector[1] = (float)normal.getY();
        nvector[2] = (float)normal.getZ();
        nvector = q.rotateVector(nvector, 0, nvector, 0);
        aV2.changeVector(nvector[0], nvector[1], nvector[2]).makeUnit().multiply(height);
        testPoint.assignTranslation(origin, aV2);
        if (Math.abs(Point.distance(patch.sphere.center, testPoint) - patch.sphere.radius) > 0.01){
            aV1.changeVector(testPoint, patch.sphere.center);
            aV1.makeUnit().multiply(sphereRadius);
            testPoint.assignTranslation(patch.sphere.center, aV1);

        }
        return testPoint;
    }

    private void closeLoop(){
        if (verbose){
            System.out.println(nextFaceID + "-th face closing loop");
        }

        if (nextFaceID >= facePool.size()){
            facePool.add(new Face(0, 0, 0));
        }
        Face nF = facePool.get(nextFaceID++);
        nF.a = e.p1._id;
        nF.b = e.p2._id;
        nF.c = e.next.p2._id;

        Surface.numoftriangles++;
        facets.remove(e);
        facets.remove(e.next);
        facets.remove(e.prev);

        nodeEdgeMap.get(e.p1.afmIdx).remove(e);
        nodeEdgeMap.get(e.p2.afmIdx).remove(e);
        nodeEdgeMap.get(e.next.p1.afmIdx).remove(e.next);
        nodeEdgeMap.get(e.next.p2.afmIdx).remove(e.next);
        nodeEdgeMap.get(e.prev.p1.afmIdx).remove(e.prev);
        nodeEdgeMap.get(e.prev.p2.afmIdx).remove(e.prev);
    }

    private void generateFaceWithNewPoint(){
        if (verbose) {
            System.out.println(nextFaceID + ". face by new face");
        }
        Point pTest = new Point(testPoint);
        pTest.afmIdx = nextNEMID;
        if (nextNEMID >= nodeEdgeMap.size()){
            nodeEdgeMap.add(new ArrayList<>());
        }
        nodeEdgeMap.get(nextNEMID).clear();
        nextNEMID++;
        nodes.add(pTest);
        pTest._id = patch.nextVertexID++;
        patch.vertices.add(pTest);
        if (nextEdgeID >= edgePool.size()){
            edgePool.add(new Edge(0, 0));
        }
        Edge leftFacet = edgePool.get(nextEdgeID++);
        if (nextEdgeID >= edgePool.size()){
            edgePool.add(new Edge(0, 0));
        }
        Edge rightFacet = edgePool.get(nextEdgeID++);
        leftFacet.v1 = e.p1.afmIdx;
        leftFacet.v2 = pTest.afmIdx;
        rightFacet.v1 = pTest.afmIdx;
        rightFacet.v2 = e.p2.afmIdx;
        leftFacet.p1 = e.p1;
        leftFacet.p2 = pTest;
        rightFacet.p1 = pTest;
        rightFacet.p2 = e.p2;
        leftFacet.prev = e.prev;
        leftFacet.prev.next = leftFacet;
        leftFacet.next = rightFacet;
        rightFacet.prev = leftFacet;
        rightFacet.next = e.next;
        rightFacet.next.prev = rightFacet;
        facets.remove(e);
        facets.add(leftFacet);
        facets.add(rightFacet);

        if (nextFaceID >= facePool.size()){
            facePool.add(new Face(0, 0, 0));
        }
        Face nF = facePool.get(nextFaceID++);
        nF.a = e.p1._id;
        nF.b = e.p2._id;
        nF.c = pTest._id;

        numOfTriangles++;
        nodeEdgeMap.get(e.p1.afmIdx).remove(e);
        nodeEdgeMap.get(e.p1.afmIdx).add(leftFacet);
        nodeEdgeMap.get(e.p2.afmIdx).remove(e);
        nodeEdgeMap.get(e.p2.afmIdx).add(rightFacet);
        nodeEdgeMap.get(pTest.afmIdx).add(leftFacet);
        nodeEdgeMap.get(pTest.afmIdx).add(rightFacet);
        e = rightFacet.next;
        leftFacet.loopID = activeLoop;
        rightFacet.loopID = activeLoop;
    }

    private void generateFaceWithPreviousEdge(){
        if (verbose) {
            System.out.println(nextFaceID + ". face constructed with prev edge");
        }
        if (nextEdgeID >= edgePool.size()){
            edgePool.add(new Edge(0, 0));
        }

        Edge newFacet = edgePool.get(nextEdgeID++);
        newFacet.v1 = e.prev.p1.afmIdx;
        newFacet.v2 = e.p2.afmIdx;
        newFacet.p1 = e.prev.p1;
        newFacet.p2 = e.p2;
        newFacet.prev = e.prev.prev;
        newFacet.prev.next = newFacet;
        newFacet.next = e.next;
        newFacet.next.prev = newFacet;
        nodeEdgeMap.get(e.prev.p1.afmIdx).remove(e.prev);
        nodeEdgeMap.get(e.prev.p1.afmIdx).add(newFacet);
        nodeEdgeMap.get(e.p1.afmIdx).remove(e.prev);
        nodeEdgeMap.get(e.p1.afmIdx).remove(e);
        nodeEdgeMap.get(e.p2.afmIdx).remove(e);
        nodeEdgeMap.get(e.p2.afmIdx).add(newFacet);

        /*
            it is safe to remove e.p1 from nodes list
         */
        nodes.remove(e.p1);
        facets.remove(e);
        facets.remove(e.prev);
        facets.add(newFacet);
        if (nextFaceID >= facePool.size()){
            facePool.add(new Face(0, 0, 0));
        }
        Face nF = facePool.get(nextFaceID++);
        nF.a = e.p1._id;
        nF.b = e.p2._id;
        nF.c = newFacet.p1._id;
        numOfTriangles++;
        e = newFacet.next;
        newFacet.loopID = activeLoop;
        faceGenerated = true;
    }

    private void generateFaceWithNextEdge(){
        if (verbose) {
            System.out.println(nextFaceID + ". face constructed with next edge");
        }
        if (nextEdgeID >= edgePool.size()){
            edgePool.add(new Edge(0, 0));
        }
        Edge newFacet = edgePool.get(nextEdgeID++);
        newFacet.v1 = e.p1.afmIdx;
        newFacet.v2 = e.next.p2.afmIdx;
        newFacet.p1 = e.p1;
        newFacet.p2 = e.next.p2;
        newFacet.prev = e.prev;
        newFacet.prev.next = newFacet;
        newFacet.next = e.next.next;
        newFacet.next.prev = newFacet;
        nodeEdgeMap.get(e.p1.afmIdx).remove(e);
        nodeEdgeMap.get(e.p1.afmIdx).add(newFacet);
        nodeEdgeMap.get(e.p2.afmIdx).remove(e);
        nodeEdgeMap.get(e.p2.afmIdx).remove(e.next);
        nodeEdgeMap.get(e.next.p2.afmIdx).remove(e.next);
        nodeEdgeMap.get(e.next.p2.afmIdx).add(newFacet);

        /*
            it is safe to remove e.p2 from nodes list
         */
        nodes.remove(e.p2);
        facets.remove(e);
        facets.remove(e.next);
        facets.add(newFacet);
        if (nextFaceID >= facePool.size()){
            facePool.add(new Face(0, 0, 0));
        }
        Face nF = facePool.get(nextFaceID++);
        nF.a = e.p1._id;
        nF.b = e.p2._id;
        nF.c = newFacet.p2._id;
        numOfTriangles++;
        e = newFacet.next;
        newFacet.loopID = activeLoop;
        faceGenerated = true;
    }

    private void generateBridgeFace(Point pt){
        if (verbose) {
            System.out.println(nextFaceID + ". face constructed with bridge edge");
        }
        List<Edge> pointEdges = nodeEdgeMap.get(pt.afmIdx);
        if (pointEdges.size() != 2){
            relevantEdges.clear();
            for (int i = 0; i < pointEdges.size(); ++i){
                Edge e = pointEdges.get(i);
                if (e.loopID == activeLoop){
                    relevantEdges.add(e);
                }
            }
            pointEdges = relevantEdges;
        }
        if (pointEdges.size() == 2) {
            Edge eNext = (pointEdges.get(0).p1 == pt) ? pointEdges.get(0) : pointEdges.get(1);
            Edge ePrev = (pointEdges.get(0).p2 == pt) ? pointEdges.get(0) : pointEdges.get(1);

            if (nextEdgeID >= edgePool.size()){
                edgePool.add(new Edge(0, 0));
            }
            Edge leftFacet = edgePool.get(nextEdgeID++);
            leftFacet.v1 = e.p1.afmIdx;
            leftFacet.v2 = eNext.p1.afmIdx;
            leftFacet.p1 = e.p1;
            leftFacet.p2 = eNext.p1;
            leftFacet.prev = e.prev;
            leftFacet.prev.next = leftFacet;
            leftFacet.next = eNext;
            leftFacet.next.prev = leftFacet;
            leftFacet.loopID = activeLoop;
            if (nextEdgeID >= edgePool.size()){
                edgePool.add(new Edge(0, 0));
            }
            Edge rightFacet = edgePool.get(nextEdgeID++);
            rightFacet.v1 = ePrev.p2.afmIdx;
            rightFacet.v2 = e.p2.afmIdx;
            rightFacet.p1 = ePrev.p2;
            rightFacet.p2 = e.p2;
            rightFacet.prev = ePrev;
            rightFacet.prev.next = rightFacet;
            rightFacet.next = e.next;
            rightFacet.next.prev = rightFacet;
            loops.add(rightFacet);
            rightFacet.loopID = numOfLoops;
            ePrev.loopID = numOfLoops;
            rightFacet.prev.loopID = numOfLoops;

            Edge aEdge = leftFacet.next;
            do {
                aEdge.loopID = activeLoop;
                aEdge = aEdge.next;
            } while (aEdge != leftFacet);

            rightFacet.loopID = numOfLoops;
            aEdge = rightFacet.next;
            do {
                aEdge.loopID = numOfLoops;
                aEdge = aEdge.next;
            } while (aEdge != rightFacet);
            numOfLoops++;

            nodeEdgeMap.get(pt.afmIdx).add(leftFacet);
            nodeEdgeMap.get(pt.afmIdx).add(rightFacet);

            nodeEdgeMap.get(e.p1.afmIdx).remove(e);
            nodeEdgeMap.get(e.p1.afmIdx).add(leftFacet);
            nodeEdgeMap.get(e.p2.afmIdx).remove(e);
            nodeEdgeMap.get(e.p2.afmIdx).add(rightFacet);
            facets.remove(e);
            facets.add(rightFacet);
            facets.add(leftFacet);
            if (nextFaceID >= facePool.size()){
                facePool.add(new Face(0, 0, 0));
            }
            Face nF = facePool.get(nextFaceID++);
            nF.a = e.p1._id;
            nF.b = e.p2._id;
            nF.c = rightFacet.p1._id;
            e = leftFacet.next;
            activeLoop = e.loopID;
            numOfTriangles++;
            faceGenerated = true;
        } else {
        }
    }

    /*
    methods used in criterions for new triangle creation
     */
    private double computeAngle(Vector v1, Vector v2, Vector normal){
        projectVectorOntoPlane(v1, normal, tau1).makeUnit();
        projectVectorOntoPlane(v2, normal, tau2).makeUnit();
        if (determinant(tau1, tau2, normal) > 0.0){
            return 2 * Math.PI - Math.acos(tau1.dotProduct(tau2));
        } else {
            return Math.acos(tau1.dotProduct(tau2));
        }
    }

    private Vector computeTriangleNormal(Point a, Point b, Point c){
        return n.assignNormalVectorOf(aV1.changeVector(b, a).makeUnit(), aV2.changeVector(c, a).makeUnit()).makeUnit();
    }

    private boolean checkForIntersectingEdges(Edge e1, Edge e2, List<Edge> toInspect, List<Edge> toIgnore){
        for (int i = 0; i < toInspect.size(); ++i){
            Edge k = toInspect.get(i);//
            if (toIgnore.contains(k)){
                continue;
            }
            if (checkForIntersectingEdges(e1, k, patch.sphere.radius, patch.sphere.center) || checkForIntersectingEdges(e2, k, patch.sphere.radius, patch.sphere.center)){
                return true;
            }
        }
        return false;
    }

    private boolean checkForIntersectingEdges(Edge e1, Edge e2, double atomRad, Point atomCenter){
        if (e1.p1 == e2.p1 || e1.p1 == e2.p2 || e1.p2 == e2.p1 || e1.p2 == e2.p2){
            return false;
        }

        _mid1.setAsMidpoint(e1.p1, e1.p2);
        _mid2.setAsMidpoint(e2.p1, e2.p2);
        double _d1 = Point.distance(e1.p1, e1.p2) * 0.5;
        double _d2 = Point.distance(e2.p1, e2.p2) * 0.5;
        _d1 += 0.1 * _d1;
        _d2 += 0.1 * _d2;

        if (Point.distance(_mid1, e2.p1) - _d1 > 0.0 && Point.distance(_mid1, e2.p2) - _d1 > 0.0
                && Point.distance(_mid2, e1.p1) - _d2 > 0.0 && Point.distance(_mid2, e1.p2) - _d2 > 0.0){
            return false;
        }

        ve1 = ve1.changeVector(e1.p1, atomCenter).makeUnit();
        ve2 = ve2.changeVector(e1 .p2, atomCenter).makeUnit();
        addv1.assignAddition(ve1, ve2).makeUnit();
        aV1.assignNormalVectorOf(ve1, ve2).makeUnit();

        v2e1 = v2e1.changeVector(e2.p1, atomCenter).makeUnit();
        v2e2 = v2e2.changeVector(e2.p2, atomCenter).makeUnit();
        addv2.assignAddition(v2e1, v2e2).makeUnit();
        if (addv1.dotProduct(addv2) < 0.0){
            return false;
        }

        if (aV1.dotProduct(v2e1) * aV1.dotProduct(v2e2) < 0){
            ro1.redefine(atomCenter, aV1);
            aV2.assignNormalVectorOf(v2e1, v2e2).makeUnit();
            ro2.redefine(atomCenter, aV2);
            if (!ro1.assignIntersectionVectorTo(intersection, ro2)){//no intersection vector was computed
                return false;
            }
            if (intersection.dotProduct(ve1) < 0 || intersection.dotProduct(ve2) < 0){
                intersection.multiply(-1);
            }
            double alpha = Math.acos(ve1.dotProduct(ve2));
            if (Math.acos(intersection.dotProduct(ve1)) - alpha < 0 && Math.acos(intersection.dotProduct(ve2)) - alpha < 0){
                return true;
            }
        }
        return false;
    }

    private double pointEdgeDistance(Point p, Edge e){
        return Point.distance(p, e.p1) + Point.distance(p, e.p2);
    }

    private Vector projectVectorOntoPlane(Vector vectorToProject, Vector normalPlaneVector, Vector tau){
        ve1.changeVector(normalPlaneVector).makeUnit();
        ve2.changeVector(vectorToProject).makeUnit();
        ve1.multiply(ve2.dotProduct(ve1));
        return tau.assignAddition(ve2, ve1.multiply(-1));
    }

    public static double determinant(Vector v1, Vector v2, Vector v3){
        return v1.getX() * v2.getY() * v3.getZ() + v1.getY() * v2.getZ() * v3.getX() + v1.getZ() * v2.getX() * v3.getY() -
                (v1.getZ() * v2.getY() * v3.getX() + v2.getZ() * v3.getY() * v1.getX() + v1.getY() * v2.getX() * v3.getZ());
    }

    private boolean hasCorrectOrientation(Point a, Point b, Point c){
        n.assignNormalVectorOf(aV1.changeVector(b, a).makeUnit(), aV2.changeVector(c, a).makeUnit()).makeUnit();
        return n.dotProduct(v3.changeVector(a, patch.sphere.center).makeUnit()) > 0.0;
    }
}

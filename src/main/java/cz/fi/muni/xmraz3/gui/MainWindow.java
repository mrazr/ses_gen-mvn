package cz.fi.muni.xmraz3.gui;


import com.jogamp.newt.MonitorDevice;
import com.jogamp.newt.Screen;
import com.jogamp.newt.event.*;
import com.jogamp.newt.opengl.GLWindow;
import com.jogamp.opengl.*;
import com.jogamp.opengl.math.Matrix4;
import com.jogamp.opengl.math.Quaternion;
import com.jogamp.opengl.math.VectorUtil;
import com.jogamp.opengl.util.Animator;
import com.jogamp.opengl.util.GLBuffers;
import com.jogamp.opengl.util.PMVMatrix;
import cz.fi.muni.xmraz3.SesConfig;
import cz.fi.muni.xmraz3.Surface;
import cz.fi.muni.xmraz3.SurfaceParser;
import cz.fi.muni.xmraz3.gui.controllers.MainPanelController;
import cz.fi.muni.xmraz3.math.Point;
import cz.fi.muni.xmraz3.math.Vector;
import cz.fi.muni.xmraz3.mesh.Boundary;
import cz.fi.muni.xmraz3.mesh.SphericalPatch;
import cz.fi.muni.xmraz3.mesh.ToroidalPatch;
import cz.fi.muni.xmraz3.utils.GLUtil;
import cz.fi.muni.xmraz3.utils.PatchUtil;
import javafx.application.Platform;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.stage.Stage;
import smile.neighbor.Neighbor;

import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

public class MainWindow implements GLEventListener, KeyListener, MouseListener{
    public static MainWindow mainWindow = null;
    public GLWindow window;
    private Animator animator;
    private int mainProgram;
    private int rectProgram;
    public Stage controlPanel;
    public BooleanProperty isPinned = new SimpleBooleanProperty(false);
    private GL4 gl;

    boolean captureMouse = false;

    //view params
    private static float[] cameraPos = new float[] {0.f, 0.f, 0.f};
    private static float[] direction = {0.f, 0.f, -1.0f};
    private float[] up = {0.0f, 1.0f, 0.0f};
    private float[] right = new float[3];
    float speed = .03f;
    float mouseSpeed = 0.02f;
    float deltaTime = 0;
    float mouseAngleX = 0.f;
    float mouseAngleY = 0.f;
    float angle = 0.0f;
    private long lastTick;
    private PMVMatrix _projMat = new PMVMatrix();
    private float[] lightPosition = new float[3];
    private float[] cameraTarget = new float[3];
    private static float[] lastCameraTarget = new float[3];
    PMVMatrix look = new PMVMatrix();

    private List<SphericalPatch> convexPatchList;
    private int hoverAtom = -1;
    public IntegerProperty selectedAtom = new SimpleIntegerProperty(1);
    public IntegerProperty selectedConcaveP = new SimpleIntegerProperty(0);
    public IntegerProperty selectedToriP = new SimpleIntegerProperty(0);
    private String strSelectedAtom = "";
    private boolean selectedExclusiveRender = false;

    private boolean drawFaces = true;
    private int drawMode = 0; // 0 -> only faces, 1 -> only lines, 2 -> faces with edges atop
    private boolean step = true;

    private List<SphericalPatch> concavePatchList;
    private boolean renderCPs = false;

    private float[] toriPatchCol = {51 / 255.f, 77 / 255.f, 177 / 255.f};
    private float[] concavePatchCol = {31 / 255.f, 143 / 255.f, 0 / 255.f};
    private float[] convexPatchCol = {197 / 255.f, 20 / 255.f, 20 / 255.f};
    private float[] selectedPatchCol = {1.f, 231 / 255.f, 76 / 255.f};
    private float[] clearColor = {1.f, 1.f, 1.f, 1.f};
    private boolean onlyCircular = false;

    private boolean mouseDown = true;
    private boolean zooming = false;
    private float zoomSpeed = 0.1f;
    private boolean cullFaces = true;
    private boolean viewPanning = false;
    private boolean mouseSelect = false;
    private int mouseSelectVbo[] = new int[1];
    private int mouseSelectVao[] = new int[1];
    private Point mouseLocation;

    static Quaternion camDir;
    static Quaternion camTar;
    private static boolean slerping = false;
    private static float slerpParam = 0.0f;

    private int[] fbo = new int[1];
    private int[] rbCol = new int[1];
    private int[] rbDep = new int[1];
    private int selectProgram = -1;
    private boolean moved = false;
    private boolean renderBuffersInit = false;
    private int uniVertexOffsets = -1;
    private int uniGlobalOffset = -1;
    private int uniStart = -1;
    private int uniEnd = -1;
    private int convexVerticesCount = 0;
    private int concaveVerticesCount = 0;
    private int toriVerticesCount = 0;
    private boolean selectInitialized = false;
    private int[] tbo = new int[1];
    private int[] tboTex = new int[1];

    private List<Integer> convexPatchesSelect = new ArrayList<>();
    private List<Integer> toriPatchesSelect = new ArrayList<>();
    private List<Integer> concavePatchesSelect = new ArrayList<>();

    private final static int CONVEX = 0;
    private final static int CONCAVE = 1;
    private final static int TORUS = 2;

    private int convexPatchesFaceCount = 0;
    private int concavePatchesFaceCount = 0;
    private int toriPatchesFaceCount = 0;

    private int[] meshVao = new int[3];
    private int[] meshVbo = new int[3];
    private int[] meshEbo = new int[2];

    private boolean convexMeshInitialized = false;
    private boolean concaveMeshInitialized = false;
    private boolean toriMeshInitialized = false;

    private int concaveFaceCountShow = 0;

    private int[] lineVao = new int[2];
    private int[] lineVbo = new int[2];
    private int[] lineEbo = new int[2];

    private int convexPatchesEdgeCount = 0;
    private int concavePatchesEdgeCount = 0;

    private AtomicBoolean convexPushData2GPU = new AtomicBoolean(false);
    private AtomicBoolean concavePushData2GPU = new AtomicBoolean(false);
    private AtomicBoolean toriPushData2GPU = new AtomicBoolean(false);

    private int[] probeVao = new int[1];
    private int[] probeVbo = new int[1];
    private int probeFaceCount = 0;
    private PMVMatrix probeScaleT = new PMVMatrix();
    private boolean renderProbe = false;
    private float probeAlpha = 0.4f;

    private List<Integer> linkNeighbors = new ArrayList<>();
    private Matrix4 selectedMeshScaleMat = new Matrix4();
    private float scaleFactor = 1.0f;
    private boolean drawModeUpdate;
    private PMVMatrix modelMAT = new PMVMatrix();
    private boolean rotating = true;
    private float ambientStrength = 0.1f;

    //pure opengl data
    private List<Integer> vaos;
    private List<Integer> vbos;
    private List<Integer> ebos;
    private int uniMeshColorLoc = -1;
    private int uniViewMatLoc = -1;
    private int uniProjMatLoc = -1;
    private int uniModelMatLoc = -1;
    private int uniAlphaLoc = -1;
    private int uniNormalColorLoc = -1;
    private int uniSelectedColorLoc = -1;
    private int uniSelectedMeshStartLoc = -1;
    private int uniSelectedMeshEndLoc = -1;
    private int uniSelectedMeshCountLoc = -1;
    private int uniAmbientStrengthLoc = -1;
    private int uniMvInverseLoc = -1;
    private int uniLightPosLoc = -1;
    private int uniCameraPosLoc = -1;
    private IntBuffer zeroes = GLBuffers.newDirectIntBuffer(new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    private AtomicBoolean stopRendering = new AtomicBoolean(true);
    private AtomicBoolean stoppedRendering = new AtomicBoolean(true);
    private boolean addToSelection;
    private boolean removeSelection;
    private boolean renderLines = false;
    private AtomicBoolean resourcesFreed = new AtomicBoolean(true);
    private int convexFaceCountShow;

    private IntBuffer mouseSelectPixelBuffer = GLBuffers.newDirectIntBuffer(1);
    private int colorChange = -1;
    private float[] newColor = new float[3];

    Quaternion axisUp = new Quaternion(0.f, 0.f, 0.f, 0.f);
    Quaternion axisRight = new Quaternion(0.f, 0.f, 0.f, 0.f);

    public void changeColor(int meshType, float r, float g, float b){
        newColor[0] = r;
        newColor[1] = g;
        newColor[2] = b;
        colorChange = meshType;
        float[] ptr = convexPatchCol;
        switch (meshType){
            case 0:
                ptr = convexPatchCol;
                break;
            case 1:
                ptr = concavePatchCol;
                break;
            case 2:
                ptr = toriPatchCol;
                break;
        }
        ptr[0] = r;
        ptr[1] = g;
        ptr[2] = b;
    }

    public void sendPatchesLists(List<SphericalPatch> convex, List<SphericalPatch> concave){
        stopRendering.set(true);
        while (!stoppedRendering.get()){
        }
        GLRunnable task = new GLRunnable() {
            @Override
            public boolean run(GLAutoDrawable glAutoDrawable) {
                freeGLResources();
                sendConvexPatchList(convex);
                sendConcavePatchList(concave);
                pushBoundariesToGPU(true);
                pushBoundariesToGPU(false);
                stopRendering.set(false);
                stoppedRendering.set(false);
                resourcesFreed.set(false);
                return true;
            }
        };
        window.invoke(false, task);
    }

    private void pushBoundariesToGPU(boolean convex){
        int totalVerticesCount = 0;
        int totalIndicesCount = 0;
        int indexOffset = 0;
        int vertexOffset = 0;
        List<Point> vertices = new ArrayList<>();
        List<Integer> indices = new ArrayList<>();
        List<SphericalPatch> patchList = (convex) ? Surface.convexPatches : Surface.triangles;
        for (SphericalPatch sp : patchList){
            sp.lineOffset = indexOffset;
            for (Boundary b : sp.boundaries){
                sp.lineCount += b.vrts.size();
                vertices.addAll(b.vrts);
                for (int i = 0; i < b.vrts.size(); ++i){
                    indices.add(i + vertexOffset);
                    if (i < b.vrts.size() - 1) {
                        indices.add(i + 1 + vertexOffset);
                    } else {
                        indices.add(vertexOffset);
                    }
                }
                vertexOffset += b.vrts.size();
            }
            indexOffset += sp.lineCount;
        }
        totalVerticesCount = vertices.size();
        totalIndicesCount = indices.size();
        if (convex){
            convexPatchesEdgeCount = indices.size() / 2;
        } else {
            concavePatchesEdgeCount = indices.size() / 2;
        }
        FloatBuffer vBuff = GLBuffers.newDirectFloatBuffer(3 * totalVerticesCount);
        IntBuffer iBuff = GLBuffers.newDirectIntBuffer(totalIndicesCount);

        vertices.forEach(p -> vBuff.put(p.getFloatData()));
        indices.forEach(i -> iBuff.put(i));

        vBuff.rewind();
        iBuff.rewind();

        gl.glBindVertexArray(lineVao[(convex) ? CONVEX : CONCAVE]);
        gl.glBindBuffer(GL4.GL_ARRAY_BUFFER, lineVbo[(convex) ? CONVEX : CONCAVE]);
        gl.glBufferData(GL4.GL_ARRAY_BUFFER, vBuff.capacity() * Float.BYTES, vBuff, GL4.GL_STATIC_DRAW);
        gl.glVertexAttribPointer(0, 3, GL4.GL_FLOAT, false, 3 * Float.BYTES, 0);
        gl.glEnableVertexAttribArray(0);
        gl.glBindBuffer(GL4.GL_ARRAY_BUFFER, 0);
        gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, lineEbo[(convex) ? CONVEX : CONCAVE]);
        gl.glBufferData(GL4.GL_ELEMENT_ARRAY_BUFFER, iBuff.capacity() * Integer.BYTES, iBuff, GL4.GL_STATIC_DRAW);
        gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
        gl.glBindVertexArray(0);
    }

    public void sendConvexPatchList(List<SphericalPatch> at){
        convexPatchList = at;
    }

    public void sendConcavePatchList(List<SphericalPatch> cp){
        concavePatchList = cp;
    }

    public void selectConvexPatchesByIDs(List<Integer> ids){
        convexPatchesSelect.clear();
        ids.forEach(i -> { if (i < convexPatchList.size()) { convexPatchesSelect.add(i);}});
    }

    public void selectConcavePatchesByIDs(List<Integer> ids){
        concavePatchesSelect.clear();
        ids.forEach(i -> { if (i < concavePatchList.size()) { concavePatchesSelect.add(i);}});
    }

    public void selectToriPatchesByIDs(List<Integer> ids){
        toriPatchesSelect.clear();
        ids.forEach(i -> { if (i < Surface.rectangles.size()) { toriPatchesSelect.add(i); }});
    }

    public void setup(){
        GLProfile profile = GLProfile.get(GLProfile.GL4);
        GLCapabilities caps = new GLCapabilities(profile);
        window = GLWindow.create(caps);
        window.setSize(800, 600);
        window.setTitle("SES");
        float[] scale = new float[2];
        window.addWindowListener(new WindowAdapter() {
            @Override
            public void windowDestroyed(WindowEvent windowEvent) {
                animator.stop();
            }
        });
        window.addWindowListener(new WindowAdapter() {
            @Override
            public void windowMoved(WindowEvent windowEvent) {
                if (controlPanel != null && isPinned.get()){
                    moveControlPanel();
                }
            }
        });
        window.addWindowListener(new WindowAdapter() {
            @Override
            public void windowResized(WindowEvent windowEvent) {
                if (controlPanel != null && isPinned.get()) {
                    moveControlPanel();
                }
            }
        });
        window.addWindowListener(new WindowAdapter() {
            @Override
            public void windowLostFocus(WindowEvent windowEvent) {
                super.windowLostFocus(windowEvent);
                window.confinePointer(false);
                window.setPointerVisible(true);
                viewPanning = false;
                zooming = false;
                captureMouse = false;
            }
        });
        isPinned.addListener(new ChangeListener<Boolean>() {
            @Override
            public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                if (newValue){
                    moveControlPanel();
                }
            }
        });
        window.addGLEventListener(this);
        window.addMouseListener(this);
        window.addKeyListener(this);
        window.setVisible(true);

        window.setResizable(true);
        Screen sc = window.getScreen();
        MonitorDevice mon = sc.getPrimaryMonitor();
        window.setPosition(mon.getViewport().getX() + mon.getViewport().getWidth() / 2 - window.getWidth() / 2,
                mon.getViewport().getY() + mon.getViewport().getHeight() / 2 - window.getHeight() / 2);

        animator = new Animator(window);
        animator.start();

        lastTick = System.currentTimeMillis();

        selectedAtom.addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                SphericalPatch sp = convexPatchList.get(newValue.intValue());
                lastCameraTarget[0] = (float)sp.sphere.center.getX();
                lastCameraTarget[1] = (float)sp.sphere.center.getY();
                lastCameraTarget[2] = (float)sp.sphere.center.getZ();
                if (!addToSelection && !removeSelection){
                    convexPatchesSelect.clear();
                    convexPatchesSelect.add(newValue.intValue());
                }
                Platform.runLater(new Runnable() {
                    @Override
                    public void run() {
                        List<String> atomInfo = new ArrayList<>();
                        SphericalPatch a = convexPatchList.get(newValue.intValue());
                        atomInfo.add("Atom radius," + Double.toString(a.sphere.radius));
                        atomInfo.add("Boundary count," + Integer.toString(a.boundaries.size()));
                        for (int i = 0; i < a.boundaries.size(); ++i){
                            atomInfo.add("Boundary " + Integer.toString(i) + "," + Integer.toString(a.boundaries.get(i).arcs.size()));
                        }
                        MainPanelController.cont.updateSelectedAtomInfo(atomInfo);
                    }
                });
            }
        });
        selectedToriP.addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if (!addToSelection && !removeSelection){
                    toriPatchesSelect.clear();
                    toriPatchesSelect.add(newValue.intValue());
                }
                ToroidalPatch tp = Surface.rectangles.get(newValue.intValue());
                lastCameraTarget[0] = (float)tp.midProbe.getX();
                lastCameraTarget[1] = (float)tp.midProbe.getY();
                lastCameraTarget[2] = (float)tp.midProbe.getZ();
            }
        });
        selectedConcaveP.addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if (!addToSelection && !removeSelection){
                    concavePatchesSelect.clear();
                    concavePatchesSelect.add(newValue.intValue());
                }
                SphericalPatch cp = concavePatchList.get(newValue.intValue());
                lastCameraTarget[0] = (float)cp.sphere.center.getX();
                lastCameraTarget[1] = (float)cp.sphere.center.getY();
                lastCameraTarget[2] = (float)cp.sphere.center.getZ();
            }
        });
        selectedMeshScaleMat.loadIdentity();
        modelMAT.glLoadIdentity();

        vaos = new ArrayList<>();
        vbos = new ArrayList<>();
        ebos = new ArrayList<>();
        MainWindow.mainWindow = this;
    }

    public void showProbe(boolean v){
        renderProbe = v;
    }

    private void moveControlPanel(){
        MainPanel.pinnedToView = false;
        controlPanel.setX(window.getX() + window.getWidth() - window.getInsets().getRightWidth());
        controlPanel.setY(window.getY() - window.getInsets().getTopHeight());
        MainPanel.pinnedToView = true;
    }

    @Override
    public void init(GLAutoDrawable glAutoDrawable) {
        gl = GLContext.getCurrentGL().getGL4();
        mainProgram = GLUtil.createShaderProgram(getClass().getClassLoader().getResourceAsStream("shaders/main.vert"), getClass().getClassLoader().getResourceAsStream("shaders/main.frag"));
        //System.out.println("Graphics vendor: " + gl.glGetString(GL.GL_VENDOR));
        uniMeshColorLoc = gl.glGetUniformLocation(mainProgram, "col");
        uniProjMatLoc = gl.glGetUniformLocation(mainProgram, "projMat");
        uniViewMatLoc = gl.glGetUniformLocation(mainProgram, "viewMat");
        uniModelMatLoc = gl.glGetUniformLocation(mainProgram, "modelMat");
        uniAlphaLoc = gl.glGetUniformLocation(mainProgram, "alpha");
        uniNormalColorLoc = gl.glGetUniformLocation(mainProgram, "normalColor");
        uniSelectedColorLoc = gl.glGetUniformLocation(mainProgram, "selectedColor");
        uniSelectedMeshStartLoc = gl.glGetUniformLocation(mainProgram, "selectedMeshStart");
        uniSelectedMeshEndLoc = gl.glGetUniformLocation(mainProgram, "selectedMeshEnd");
        uniAmbientStrengthLoc = gl.glGetUniformLocation(mainProgram, "ambientStrength");
        uniSelectedMeshCountLoc = gl.glGetUniformLocation(mainProgram, "selectedMeshCount");
        uniMvInverseLoc = gl.glGetUniformLocation(mainProgram, "mvInverse");
        uniCameraPosLoc = gl.glGetUniformLocation(mainProgram, "cameraPos");
        uniLightPosLoc = gl.glGetUniformLocation(mainProgram, "lightPos");
        selectProgram = GLUtil.createShaderProgram(getClass().getClassLoader().getResourceAsStream("shaders/selectShader.vert"), getClass().getClassLoader().getResourceAsStream("shaders/selectShader.frag"));
        uniEnd = gl.glGetUniformLocation(selectProgram, "end");
        uniStart = gl.glGetUniformLocation(selectProgram, "start");
        uniGlobalOffset = gl.glGetUniformLocation(selectProgram, "globalOffset");
        uniVertexOffsets = gl.glGetUniformLocation(selectProgram, "u_offset_Tex");
        right = VectorUtil.crossVec3(right, direction, up);
        up = VectorUtil.normalizeVec3(up);
        right = VectorUtil.normalizeVec3(right);
        direction = VectorUtil.normalizeVec3(direction);
        gl.glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
        gl.glEnable(GL.GL_DEPTH_TEST);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glEnable(GL.GL_CULL_FACE);
        gl.glLineWidth(0.5f);
        gl.glCullFace(GL.GL_BACK);
        float aspect = 800 / (float)600;
        _projMat.glLoadIdentity();
        _projMat.gluPerspective(60.f, aspect, 0.1f, 1000.f);
        lastTick = System.currentTimeMillis();
        mouseLocation = new Point(0, 0, 0);
        gl.glGenVertexArrays(1, mouseSelectVao, 0);
        gl.glGenBuffers(1, mouseSelectVbo, 0);
        FloatBuffer bff = GLBuffers.newDirectFloatBuffer(6);
        bff.put(new float[]{0.f, 0.f, 0.f});
        bff.rewind();
        gl.glBindVertexArray(mouseSelectVao[0]);
        gl.glBindBuffer(GL.GL_ARRAY_BUFFER, mouseSelectVbo[0]);
        gl.glBufferData(GL.GL_ARRAY_BUFFER, bff.capacity() * Float.BYTES, bff, GL.GL_DYNAMIC_DRAW);
        gl.glVertexAttribPointer(0, 3, GL.GL_FLOAT, false, 3 * Float.BYTES, 0);
        gl.glEnableVertexAttribArray(0);
        gl.glBindBuffer(GL.GL_ARRAY_BUFFER, 0);
        gl.glBindVertexArray(0);

        gl.glGenFramebuffers(1, fbo, 0);
        gl.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo[0]);
        gl.glGenRenderbuffers(1, rbCol, 0);
        gl.glGenRenderbuffers(1, rbDep, 0);
        gl.glBindRenderbuffer(GL.GL_RENDERBUFFER, rbCol[0]);
        gl.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL4.GL_R32I, window.getWidth(), window.getHeight());
        gl.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL4.GL_COLOR_ATTACHMENT0, GL4.GL_RENDERBUFFER, rbCol[0]);
        gl.glBindRenderbuffer(GL.GL_RENDERBUFFER, rbDep[0]);
        gl.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL4.GL_DEPTH_COMPONENT24, window.getWidth(), window.getHeight());
        gl.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL4.GL_DEPTH_ATTACHMENT, GL4.GL_RENDERBUFFER, rbDep[0]);
        gl.glBindRenderbuffer(GL.GL_RENDERBUFFER, 0);
        gl.glBindFramebuffer(GL4.GL_FRAMEBUFFER, 0);
        int[] bobjects = GLUtil.loadSphere(getClass().getClassLoader().getResourceAsStream("misc/sphere3.obj"), gl);
        this.probeVao[0] = bobjects[0];
        this.probeVbo[0] = bobjects[1];
        this.probeFaceCount = bobjects[2];
        this.probeScaleT.glLoadIdentity();
        float r = (float)Double.longBitsToDouble(Surface.probeRadius.get());
        this.probeScaleT.glPushMatrix();
        camDir = new Quaternion(direction[0] ,direction[1], direction[2], 0.f);
        camDir.normalize();
        camTar = new Quaternion((float) Surface.centerOfgravity.getX() - cameraPos[0], (float) Surface.centerOfgravity.getY() - cameraPos[1], (float) Surface.centerOfgravity.getZ() - cameraPos[2], 0.f);
        camTar.normalize();
        slerping = true;
        slerpParam = 0;
        cameraTarget[0] = 0.f;
        cameraTarget[1] = 0.f;
        cameraTarget[2] = -2.f;
        lastCameraTarget[0] = cameraTarget[0];
        lastCameraTarget[1] = cameraTarget[1];
        lastCameraTarget[2] = cameraTarget[2];
        lightPosition[0] = cameraPos[0];
        lightPosition[1] = cameraPos[1];
        lightPosition[2] = cameraPos[2];
        right = new float[]{1.0f, 0.f, 0.f};
        up = VectorUtil.crossVec3(up, right, direction);
    }

    @Override
    public void dispose(GLAutoDrawable glAutoDrawable) {
        stoppedRendering.set(true);
        freeGLResources();
        animator.stop();
    }

    private void deleteRenderBuffers(){
        gl.glDeleteRenderbuffers(1, rbCol, 0);
        gl.glDeleteRenderbuffers(1, rbDep, 0);
    }

    private void constructNewRenderBuffers(){
        if (renderBuffersInit){
            deleteRenderBuffers();
        }
        gl.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo[0]);
        gl.glGenRenderbuffers(1, rbCol, 0);
        gl.glGenRenderbuffers(1, rbDep, 0);
        gl.glBindRenderbuffer(GL.GL_RENDERBUFFER, rbCol[0]);
        gl.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL4.GL_R32I, window.getWidth(), window.getHeight());
        gl.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL4.GL_COLOR_ATTACHMENT0, GL4.GL_RENDERBUFFER, rbCol[0]);
        gl.glBindRenderbuffer(GL.GL_RENDERBUFFER, rbDep[0]);
        gl.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL4.GL_DEPTH_COMPONENT24, window.getWidth(), window.getHeight());
        gl.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL4.GL_DEPTH_ATTACHMENT, GL4.GL_RENDERBUFFER, rbDep[0]);

        if (gl.glCheckFramebufferStatus(GL4.GL_FRAMEBUFFER) == GL4.GL_FRAMEBUFFER_COMPLETE){
            renderBuffersInit = true;
        } else {
            renderBuffersInit = false;
            deleteRenderBuffers();
        }
        gl.glBindRenderbuffer(GL.GL_RENDERBUFFER, 0);
        gl.glBindFramebuffer(GL4.GL_FRAMEBUFFER, 0);

    }
    private void renderProbeAndNeighborProbes(){
        if (renderProbe) {
            gl.glUseProgram(mainProgram);
            gl.glEnable(GL.GL_BLEND);
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glUniform3f(uniMeshColorLoc, 1.f, 1.f, 1.f);
            gl.glUniform1f(uniAlphaLoc, probeAlpha);
            if (linkNeighbors.size() == 0) {
                float r = (float)SesConfig.probeRadius;
                for (int j = 0; j < concavePatchesSelect.size(); ++j){
                    Integer i = concavePatchesSelect.get(j);
                    SphericalPatch cp2 = concavePatchList.get(i);
                    float[] c2 = cp2.sphere.center.getFloatData();
                    probeScaleT.glPushMatrix();
                    probeScaleT.glTranslatef(c2[0], c2[1], c2[2]);
                    probeScaleT.glScalef(r, r, r);
                    gl.glUniformMatrix4fv(uniModelMatLoc, 1, false, probeScaleT.glGetMatrixf());
                    gl.glBindVertexArray(probeVao[0]);
                    gl.glDrawArrays(GL4.GL_TRIANGLES, 0, 3 * probeFaceCount);
                    gl.glBindVertexArray(0);
                    probeScaleT.glPopMatrix();
                }

            } else {
                for (int j = 0; j < linkNeighbors.size(); ++j) {
                    Integer i = linkNeighbors.get(j);
                    SphericalPatch cp = concavePatchList.get(i);
                    probeScaleT.glPushMatrix();
                    float r = (float) Double.longBitsToDouble(Surface.probeRadius.get());

                    float[] center = cp.sphere.center.getFloatData();
                    probeScaleT.glTranslatef(center[0], center[1], center[2]);
                    probeScaleT.glScalef(r, r, r);
                    gl.glUniformMatrix4fv(uniModelMatLoc, 1, false, probeScaleT.glGetMatrixf());
                    gl.glBindVertexArray(probeVao[0]);
                    gl.glDrawArrays(GL.GL_TRIANGLES, 0, 3 * probeFaceCount);
                    gl.glBindVertexArray(0);
                    probeScaleT.glPopMatrix();
                }
            }
            gl.glDisable(GL.GL_BLEND);
        }
    }
    @Override
    public void display(GLAutoDrawable glAutoDrawable) {
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
        updateCamera();
        if (!stopRendering.get()){
            if (mouseSelect && moved){
                moved = false;
                gl.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo[0]);
                mouseSelectPixelBuffer.clear();
                gl.glReadPixels((int)mouseLocation.x, window.getHeight() - 1 - (int)mouseLocation.y, 1, 1, GL4.GL_RED_INTEGER, GL4.GL_INT, mouseSelectPixelBuffer);
                hoverAtom = mouseSelectPixelBuffer.get();
                gl.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0);
            }
            gl.glUseProgram(mainProgram);
            int lightColor_loc = gl.glGetUniformLocation(mainProgram, "lightColor");
            gl.glUniformMatrix4fv(uniProjMatLoc, 1, false, _projMat.glGetMatrixf());
            gl.glUniformMatrix4fv(uniViewMatLoc, 1, false, look.glGetMatrixf());
            gl.glUniformMatrix4fv(uniModelMatLoc, 1, false, modelMAT.glGetMatrixf());
            gl.glUniformMatrix4fv(uniMvInverseLoc, 1, false, modelMAT.glGetMvitMatrixf());
            gl.glUniform3f(lightColor_loc, 1.f, 1.f, 1.f);
            gl.glUniform3fv(uniLightPosLoc, 1, lightPosition, 0);
            gl.glUniform1f(uniAlphaLoc, 1.f);
            drawConvex();
            drawConcave();
            drawTori();
            renderProbeAndNeighborProbes();
        } else {
            stoppedRendering.set(true);
        }
        gl.glUseProgram(0);
        gl.glBindVertexArray(0);
        deltaTime = (float)(System.currentTimeMillis() - lastTick);
        angle += deltaTime * 0.5f;
        lastTick = System.currentTimeMillis();
    }

    @Override
    public void reshape(GLAutoDrawable glAutoDrawable, int i, int i1, int i2, int i3) {
        float aspect = (float)i2 / (float)i3;
        _projMat.glLoadIdentity();
        _projMat.gluPerspective(60.f, aspect, 0.1f, 1000.f);
        constructNewRenderBuffers();
        updateCamera();
    }

    private int cforward = 0;
    private int cright = 0;

    private void updateCamera(){
        cameraPos[0] += direction[0] * cforward * speed * deltaTime;
        cameraPos[1] += direction[1] * cforward * speed * deltaTime;
        cameraPos[2] += direction[2] * cforward * speed * deltaTime;

        cameraPos[0] += right[0] * cright * speed * deltaTime;
        cameraPos[1] += right[1] * cright * speed * deltaTime;
        cameraPos[2] += right[2] * cright * speed * deltaTime;

        lightPosition[0] = cameraPos[0];
        lightPosition[1] = cameraPos[1];
        lightPosition[2] = cameraPos[2];
        if (slerping && slerpParam < 1.f){
            if (camDir.dot(camTar) < 0.0f){
                camDir.scale(-1.f);
            }
            Quaternion slerpq = new Quaternion(0.f, 0.f, 0.f, 0.f);
            slerpq.setSlerp(camDir, camTar, slerpParam);
            slerpParam += 0.01f * deltaTime;
            direction[0] = slerpq.getX();
            direction[1] = slerpq.getY();
            direction[2] = slerpq.getZ();
            direction = VectorUtil.normalizeVec3(direction);
            right = VectorUtil.crossVec3(right, new float[]{0.f, 1.f, 0.f}, direction);
            up = VectorUtil.crossVec3(up, right, direction);
        } else {
            slerping = false;
            slerpParam = 0.f;
        }

        look.glLoadIdentity();
        look.gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2],
                cameraPos[0] + direction[0], cameraPos[1] + direction[1], + cameraPos[2] + direction[2],
                up[0], up[1], up[2]);

        if (!stopRendering.get() && mouseSelect){
            if (convexPatchesFaceCount == 0 || concavePatchesFaceCount == 0 || toriPatchesFaceCount == 0){
                return;
            }
            gl.glUseProgram(selectProgram);
            if (!selectInitialized){
                IntBuffer boffsets = GLBuffers.newDirectIntBuffer(convexPatchList.size() + concavePatchList.size() + Surface.rectangles.size());
                int accumulator = 0;
                for (SphericalPatch a : convexPatchList){
                    boffsets.put(accumulator);
                    accumulator += a.vertices.size();
                }
                convexVerticesCount = accumulator;
                for (SphericalPatch cp : concavePatchList){
                    boffsets.put(accumulator);
                    accumulator += cp.vertices.size();
                    concaveVerticesCount += cp.vertices.size();
                }
                for (ToroidalPatch tp : Surface.rectangles){
                    if (tp.faces == null){
                        continue;
                    }
                    boffsets.put(accumulator);
                    accumulator += tp.faces.length;// * 3;//tp.vrts.size() / 2;
                    toriVerticesCount += tp.faces.length;//tp.vrts.size() / 2;
                }
                //toriVerticesCount = accumulator;
                boffsets.rewind();
                gl.glGenBuffers(1, tbo, 0);
                gl.glBindBuffer(GL4.GL_TEXTURE_BUFFER, tbo[0]);
                gl.glBufferData(GL4.GL_TEXTURE_BUFFER, boffsets.capacity() * Integer.BYTES, boffsets, GL4.GL_STATIC_DRAW);
                gl.glGenTextures(1, tboTex, 0);
                gl.glBindBuffer(GL4.GL_TEXTURE_BUFFER, 0);
                selectInitialized = true;
            }
            gl.glUseProgram(selectProgram);
            gl.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo[0]);
            gl.glEnable(GL.GL_DEPTH_TEST);
            gl.glClearColor(-1.f, -1.f, -1.f, -1.f);
            gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);

            gl.glActiveTexture(GL4.GL_TEXTURE0);
            gl.glBindTexture(GL4.GL_TEXTURE_BUFFER, tboTex[0]);
            gl.glTexBuffer(GL4.GL_TEXTURE_BUFFER, GL4.GL_R32I, tbo[0]);
            gl.glUniform1i(uniVertexOffsets, 0);


            int viewLoc = gl.glGetUniformLocation(selectProgram, "viewMat");
            int projLoc = gl.glGetUniformLocation(selectProgram, "proj_matrix");
            int modelLoc = gl.glGetUniformLocation(selectProgram, "modelMat");

            gl.glUniformMatrix4fv(projLoc, 1, false, _projMat.glGetMatrixf());
            gl.glUniformMatrix4fv(viewLoc, 1, false, look.glGetMatrixf());
            gl.glUniformMatrix4fv(modelLoc, 1, false, modelMAT.glGetMatrixf());

            gl.glUniform1i(uniStart, 0);
            gl.glUniform1i(uniEnd, convexPatchList.size());
            gl.glUniform1i(uniGlobalOffset, 0);
            gl.glBindVertexArray(meshVao[CONVEX]);
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[CONVEX]);
            gl.glDrawElements(GL4.GL_TRIANGLES, 3 * convexPatchesFaceCount, GL4.GL_UNSIGNED_INT, 0);

            gl.glFrontFace(GL4.GL_CW);
            gl.glUniform1i(uniGlobalOffset, convexVerticesCount);
            gl.glUniform1i(uniStart, convexPatchList.size());
            gl.glUniform1i(uniEnd, convexPatchList.size() + concavePatchList.size());
            gl.glBindVertexArray(meshVao[CONCAVE]);
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[CONCAVE]);
            gl.glDrawElements(GL4.GL_TRIANGLES, 3 * concavePatchesFaceCount, GL4.GL_UNSIGNED_INT, 0);
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
            gl.glFrontFace(GL4.GL_CCW);

            gl.glUniform1i(uniStart, concavePatchList.size() + convexPatchList.size());
            gl.glUniform1i(uniEnd, concavePatchList.size() + convexPatchList.size() + Surface.rectangles.size());
            gl.glUniform1i(uniGlobalOffset, convexVerticesCount + concaveVerticesCount);
            gl.glBindVertexArray(meshVao[TORUS]);
            gl.glDrawArrays(GL4.GL_TRIANGLES, 0, 3 * toriPatchesFaceCount);
            gl.glBindVertexArray(0);
            gl.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0);
            gl.glUseProgram(mainProgram);
            gl.glClearColor(1.f, 1.f, 1.f, 1.f);
        }
    }

    private IntBuffer selectStart = GLBuffers.newDirectIntBuffer(20);
    private IntBuffer selectEnd = GLBuffers.newDirectIntBuffer(20);

    private void drawConvex(){
        gl.glUniform3fv(uniNormalColorLoc, 1, convexPatchCol, 0);
        gl.glUniform3fv(uniSelectedColorLoc, 1, selectedPatchCol, 0);
        gl.glUniform1f(uniAmbientStrengthLoc, ambientStrength);
        if (convexMeshInitialized) {
            if (selectedExclusiveRender){
                if (convexPatchesSelect.size() > 0) {
                    gl.glUniform1i(uniSelectedMeshCountLoc, 0);
                    gl.glUniform3fv(uniNormalColorLoc, 1, convexPatchCol, 0);
                    gl.glBindVertexArray(meshVao[CONVEX]);
                    gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[CONVEX]);
                    for (int i = 0; i < convexPatchesSelect.size(); ++i) {
                        Integer j = convexPatchesSelect.get(i);
                        SphericalPatch a = convexPatchList.get(j);
                        int off = (convexFaceCountShow > a.faces.length / 3) ? a.faces.length : convexFaceCountShow * 3;
                        gl.glDrawElements(GL4.GL_TRIANGLES,  a.faces.length - off, GL4.GL_UNSIGNED_INT, a.eboOffset * Integer.BYTES);
                    }
                }
            } else {
                gl.glUniform1iv(uniSelectedMeshStartLoc, 20, zeroes);
                gl.glUniform1iv(uniSelectedMeshEndLoc, 20, zeroes);
                selectStart.clear();
                selectEnd.clear();
                for (int i = 0; i < convexPatchesSelect.size(); ++i) {
                    Integer j = convexPatchesSelect.get(i);
                    SphericalPatch a = convexPatchList.get(j);
                    selectStart.put(a.vboOffset);
                    selectEnd.put(a.vboOffset + a.vertices.size());
                }
                if (convexPatchesSelect.size() == 0) {
                    for (int i = 0; i < 20; ++i) {
                        selectStart.put(-1);
                    }
                }
                selectStart.rewind();
                selectEnd.rewind();
                int buffCap = convexPatchesSelect.size();
                gl.glUniform1i(uniSelectedMeshCountLoc, buffCap); //convexPatchesSelect.size());
                gl.glUniform1iv(uniSelectedMeshStartLoc, buffCap, selectStart);
                gl.glUniform1iv(uniSelectedMeshEndLoc, buffCap, selectEnd);
                gl.glUniform3fv(uniNormalColorLoc, 1, convexPatchCol, 0);
                gl.glBindVertexArray(meshVao[CONVEX]);
                gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[CONVEX]);
                gl.glDrawElements(GL4.GL_TRIANGLES, 3 * convexPatchesFaceCount, GL4.GL_UNSIGNED_INT, 0);
                gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
                gl.glBindVertexArray(0);
            }
        }
        if (renderLines) {
            gl.glUniform1f(uniAmbientStrengthLoc, 1.f);
            gl.glBindVertexArray(lineVao[CONVEX]);
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, lineEbo[CONVEX]);
            if (selectedExclusiveRender){
                for (int i = 0; i < convexPatchesSelect.size(); ++i){
                    Integer j = convexPatchesSelect.get(i);
                    SphericalPatch a = convexPatchList.get(j);
                    gl.glDrawElements(GL4.GL_LINES, 2 * a.lineCount, GL4.GL_UNSIGNED_INT, 2 * a.lineOffset * Integer.BYTES);
                }
            } else {
                gl.glDrawElements(GL4.GL_LINES, 2 * convexPatchesEdgeCount, GL4.GL_UNSIGNED_INT, 0);
            }
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
            gl.glBindVertexArray(0);
        }
    }

    private void drawConcave(){
        gl.glUniform3fv(uniNormalColorLoc, 1, concavePatchCol, 0);
        gl.glUniform3fv(uniSelectedColorLoc, 1, selectedPatchCol, 0);
        gl.glUniform1f(uniAmbientStrengthLoc, ambientStrength);
        if (concaveMeshInitialized) {
            if (selectedExclusiveRender){
                if (concavePatchesSelect.size() > 0) {
                    gl.glFrontFace(GL4.GL_CW);
                    gl.glUniform1i(uniSelectedMeshCountLoc, 0);
                    gl.glUniform3fv(uniNormalColorLoc, 1, concavePatchCol, 0);
                    gl.glBindVertexArray(meshVao[CONCAVE]);
                    gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[CONCAVE]);
                    for (int i = 0; i < concavePatchesSelect.size(); ++i) {
                        Integer j = concavePatchesSelect.get(i);
                        SphericalPatch a = concavePatchList.get(j);
                        int off = (concaveFaceCountShow > a.faces.length / 3) ? a.faces.length / 3 : concaveFaceCountShow * 3;
                        gl.glDrawElements(GL4.GL_TRIANGLES, a.faces.length - off, GL4.GL_UNSIGNED_INT, a.eboOffset * Integer.BYTES);
                    }
                    gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
                    gl.glBindVertexArray(0);
                    gl.glFrontFace(GL4.GL_CCW);
                }
            } else {
                gl.glUniform1iv(uniSelectedMeshStartLoc, 20, zeroes);
                gl.glUniform1iv(uniSelectedMeshEndLoc, 20, zeroes);
                selectStart.clear();
                selectEnd.clear();
                for (int i = 0; i < concavePatchesSelect.size(); ++i) {
                    Integer j = concavePatchesSelect.get(i);
                    SphericalPatch a = concavePatchList.get(j);
                    selectStart.put(a.vboOffset);
                    selectEnd.put(a.vboOffset + a.vertices.size());
                }
                if (concavePatchesSelect.size() == 0) {
                    for (int i = 0; i < 20; ++i) {
                        selectStart.put(-1);
                    }
                }
                selectStart.rewind();
                selectEnd.rewind();
                int  buffCap = concavePatchesSelect.size();
                gl.glFrontFace(GL4.GL_CW);
                gl.glUniform1i(uniSelectedMeshCountLoc, buffCap); //concavePatchesSelect.size());
                gl.glUniform1iv(uniSelectedMeshStartLoc, buffCap, selectStart);
                gl.glUniform1iv(uniSelectedMeshEndLoc, buffCap, selectEnd);
                gl.glFrontFace(GL4.GL_CW);
                gl.glUniform3fv(uniNormalColorLoc, 1, concavePatchCol, 0);
                gl.glBindVertexArray(meshVao[CONCAVE]);
                gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[CONCAVE]);
                gl.glDrawElements(GL4.GL_TRIANGLES, 3 * concavePatchesFaceCount, GL4.GL_UNSIGNED_INT, 0);
                gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
                gl.glBindVertexArray(0);
                gl.glFrontFace(GL4.GL_CCW);
            }
        }
        if (renderLines) {
            gl.glUniform1f(uniAmbientStrengthLoc, .5f);
            gl.glBindVertexArray(lineVao[CONCAVE]);
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, lineEbo[CONCAVE]);
            if (selectedExclusiveRender){
                for (int i = 0; i < concavePatchesSelect.size(); ++i){
                    Integer j = concavePatchesSelect.get(i);
                    SphericalPatch cp = concavePatchList.get(j);
                    gl.glDrawElements(GL4.GL_LINES, 2 * cp.lineCount, GL4.GL_UNSIGNED_INT, 2 * cp.lineOffset * Integer.BYTES);
                }
            } else {
                gl.glDrawElements(GL4.GL_LINES, 2 * concavePatchesEdgeCount, GL4.GL_UNSIGNED_INT, 0);
            }
            gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
            gl.glBindVertexArray(0);
        }
    }

    private void drawTori(){
        if (toriMeshInitialized) {
            if (selectedExclusiveRender){
                if (toriPatchesSelect.size() == 0){
                    return;
                }
                gl.glUniform1i(uniSelectedMeshCountLoc, 0);
                gl.glUniform3fv(uniNormalColorLoc, 1, toriPatchCol, 0);
                gl.glBindVertexArray(meshVao[TORUS]);
                for (int i = 0; i < toriPatchesSelect.size(); ++i){
                    Integer j = toriPatchesSelect.get(i);
                    ToroidalPatch rp = Surface.rectangles.get(j);
                    gl.glDrawArrays(GL4.GL_TRIANGLES, rp.vboOffset, rp.faces.length);
                }
                gl.glBindVertexArray(0);
                return;
            }
            gl.glUniform1iv(uniSelectedMeshStartLoc, 20, zeroes);
            gl.glUniform1iv(uniSelectedMeshEndLoc, 20, zeroes);
            selectStart.clear();
            selectEnd.clear();
            for (int i = 0; i < toriPatchesSelect.size(); ++i) {
                Integer j = toriPatchesSelect.get(i);
                ToroidalPatch a = Surface.rectangles.get(j);
                selectStart.put(a.vboOffset);
                selectEnd.put(a.vboOffset + (a.faces.length));//a.vrts.size() / 2));
            }
            if (toriPatchesSelect.size() == 0) {
                for (int i = 0; i < 20; ++i) {
                    selectStart.put(-1);
                }
            }
            selectStart.rewind();
            selectEnd.rewind();
            int buffCap = toriPatchesSelect.size();//(toriPatchesSelect.size() > 0) ? toriPatchesSelect.size() : 20;
            gl.glUniform1f(uniAmbientStrengthLoc, ambientStrength);
            gl.glUniform1i(uniSelectedMeshCountLoc, buffCap);//toriPatchesSelect.size());
            gl.glUniform1iv(uniSelectedMeshStartLoc, buffCap, selectStart);
            gl.glUniform1iv(uniSelectedMeshEndLoc, buffCap, selectEnd);
            gl.glUniform3fv(uniNormalColorLoc, 1, toriPatchCol, 0);
            gl.glBindVertexArray(meshVao[TORUS]);
            gl.glDrawArrays(GL4.GL_TRIANGLES, 0, 3 * toriPatchesFaceCount);
            gl.glBindVertexArray(0);
        }
    }

    private void pushMeshToGPU(List<Point> vrtsNnormals, List<Integer> indices, int bufferObjectsIdx){
        stopRendering.set(true);
        FloatBuffer VNBuffer = GLBuffers.newDirectFloatBuffer(vrtsNnormals.size() * 3);
        for (Point p : vrtsNnormals){
            VNBuffer.put(p.getFloatData());
        }
        IntBuffer indBuffer = GLBuffers.newDirectIntBuffer(indices.size());
        for (Integer i : indices){
            indBuffer.put(i);
        }
        VNBuffer.rewind();
        indBuffer.rewind();
        gl.glBindVertexArray(meshVao[bufferObjectsIdx]);
        gl.glBindBuffer(GL4.GL_ARRAY_BUFFER, meshVbo[bufferObjectsIdx]);
        gl.glBufferData(GL4.GL_ARRAY_BUFFER, VNBuffer.capacity() * Float.BYTES, VNBuffer, GL4.GL_STATIC_DRAW);
        gl.glVertexAttribPointer(0, 3, GL4.GL_FLOAT, false, 6 * Float.BYTES, 0);
        gl.glEnableVertexAttribArray(0);
        gl.glVertexAttribPointer(1, 3, GL4.GL_FLOAT, false, 6 * Float.BYTES, 3 * Float.BYTES);
        gl.glEnableVertexAttribArray(1);
        gl.glBindBuffer(GL4.GL_ARRAY_BUFFER, 0);
        gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, meshEbo[bufferObjectsIdx]);
        gl.glBufferData(GL4.GL_ELEMENT_ARRAY_BUFFER, indBuffer.capacity() * Integer.BYTES, indBuffer, GL4.GL_STATIC_DRAW);
        gl.glBindBuffer(GL4.GL_ELEMENT_ARRAY_BUFFER, 0);
        gl.glBindVertexArray(0);
        stopRendering.set(false);
        if (bufferObjectsIdx == CONCAVE){
            concaveMeshInitialized = true;
        } else {
            convexMeshInitialized = true;
        }
    }

    private void pushConvexPatchesToGPU(){
        List<Point> vrtsNnormals = new ArrayList<>();
        List<Integer> indices = new ArrayList<>();
        int vboOffset = 0;
        int eboOffset = 0;
        int faceCount = 0;
        for (SphericalPatch a : Surface.convexPatches){
            if (!a.meshed){
                continue;
            }
            for (Point v : a.vertices){
                Point n = new Point(Point.subtractPoints(v, a.sphere.center).makeUnit().getFloatData());
                vrtsNnormals.add(v);
                vrtsNnormals.add(n);
            }
            for (int i = 0; i < a.faces.length; ++i){
                indices.add(a.faces[i] + vboOffset);
            }
            a.vboOffset = vboOffset;
            a.eboOffset = eboOffset;
            vboOffset += a.vertices.size();
            eboOffset += a.faces.length;
            faceCount += a.faces.length / 3;
        }
        convexPatchesFaceCount = indices.size() / 3;
        pushMeshToGPU(vrtsNnormals, indices, CONVEX);
        convexPushData2GPU.set(false);
        if (SesConfig.verbose) {
            System.out.println("Number of triangles on convex patches: " + faceCount);
        }
}

    private void pushConcavePatchesToGPU(){
        List<Point> vrtsNnormals = new ArrayList<>();
        List<Integer> indices = new ArrayList<>();
        int vboOffset = 0;
        int eboOffset = 0;
        int faceCount = 0;
        for (SphericalPatch cp : Surface.triangles){
            if (!cp.meshed){
                continue;
            }
            for (Point v : cp.vertices){
                Point n = new Point(Point.subtractPoints(cp.sphere.center, v).makeUnit().getFloatData());
                vrtsNnormals.add(v);
                vrtsNnormals.add(n);
            }
            for (int i = 0; i < cp.faces.length; i++){ //Face f : cp.faces){
                indices.add(cp.faces[i] + vboOffset);
            }
            cp.vboOffset = vboOffset;
            cp.eboOffset = eboOffset;
            vboOffset += cp.vertices.size();
            eboOffset += cp.faces.length;
            faceCount += cp.faces.length / 3;
        }
        concavePatchesFaceCount = indices.size() / 3;
        pushMeshToGPU(vrtsNnormals, indices, CONCAVE);
        concavePushData2GPU.set(false);
        if (SesConfig.verbose) {
            System.out.println("Number of triangles on concave patches: " + faceCount);
        }
    }

    public void pushConvex(){
        GLRunnable r = new GLRunnable() {
            @Override
            public boolean run(GLAutoDrawable glAutoDrawable) {
                pushConvexPatchesToGPU();
                return true;
            }
        };
        window.invoke(true, r);
    }

    public void pushConcave(){
        GLRunnable r = new GLRunnable() {
            @Override
            public boolean run(GLAutoDrawable glAutoDrawable) {
                pushConcavePatchesToGPU();
                return true;
            }
        };
        window.invoke(true, r);
    }


    private void _pushToriMeshToGPU(){
        stopRendering.set(true);
        FloatBuffer _vrtsNormals = GLBuffers.newDirectFloatBuffer(Surface.toriFacesCount * 3 * 2 * 3);
        int vboOffset = 0;
        int faceCount = 0;
        int _tmpOffset = 0;
        Point _p;
        Point _probe;
        Vector _n = new Vector(0, 0, 0);
        for (ToroidalPatch tp : Surface.rectangles){
            if (!tp.valid){
                continue;
            }
            _tmpOffset = 0;
            for (int i = 0; i < tp.faces.length; i ++){//Face f : tp.faces){
               _p = tp.vertices.get(tp.faces[i]);
               _probe = tp.probes[PatchUtil.getTorusProbeIdx(tp, tp.faces[i])];
               _vrtsNormals.put((float)_p.getX());
               _vrtsNormals.put((float)_p.getY());
               _vrtsNormals.put((float)_p.getZ());

               _n.changeVector(_probe, _p).makeUnit();
               _vrtsNormals.put((float)_n.getX());
               _vrtsNormals.put((float)_n.getY());
               _vrtsNormals.put((float)_n.getZ());
               _tmpOffset += 1;
            }
            tp.vboOffset = vboOffset;
            vboOffset += _tmpOffset;
            faceCount += tp.faces.length / 3;
        }
        _vrtsNormals.rewind();
        toriPatchesFaceCount = _vrtsNormals.capacity() / 18;
        gl.glBindVertexArray(meshVao[TORUS]);
        gl.glBindBuffer(GL4.GL_ARRAY_BUFFER, meshVbo[TORUS]);
        gl.glBufferData(GL4.GL_ARRAY_BUFFER, _vrtsNormals.capacity() * Float.BYTES, _vrtsNormals, GL4.GL_STATIC_DRAW);
        gl.glVertexAttribPointer(0, 3, GL4.GL_FLOAT, false, 6 * Float.BYTES, 0);
        gl.glEnableVertexAttribArray(0);
        gl.glVertexAttribPointer(1, 3, GL4.GL_FLOAT, false, 6 * Float.BYTES, 3 * Float.BYTES);
        gl.glEnableVertexAttribArray(1);
        gl.glBindBuffer(GL4.GL_ARRAY_BUFFER, 0);
        gl.glBindVertexArray(0);
        stopRendering.set(false);
        toriPushData2GPU.set(false);
        toriMeshInitialized = true;
        if (SesConfig.verbose) {
            System.out.println("Number of triangles toroidal patches: " + faceCount);
        }
        _vrtsNormals.clear();
    }

    public void pushTori(){
        GLRunnable r = new GLRunnable() {
            @Override
            public boolean run(GLAutoDrawable glAutoDrawable) {
                _pushToriMeshToGPU();
                return true;
            }
        };
        window.invoke(true, r);
    }

    @Override
    public void keyPressed(KeyEvent keyEvent) {
        if (keyEvent.getKeyCode() == KeyEvent.VK_F1){
            concaveFaceCountShow = (concaveFaceCountShow > 0) ? concaveFaceCountShow - 1 : concaveFaceCountShow;
            convexFaceCountShow = (convexFaceCountShow > 0) ? convexFaceCountShow - 1 : convexFaceCountShow;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_F2){
            concaveFaceCountShow++;
            convexFaceCountShow++;
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_F3){
            concaveFaceCountShow = 0;
            convexFaceCountShow = 0;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_F4){
            concaveFaceCountShow = (concavePatchesSelect.size() > 0) ? Surface.triangles.get(concavePatchesSelect.get(0)).faces.length / 3: 0;
            convexFaceCountShow = (convexPatchesSelect.size() > 0) ? Surface.convexPatches.get(convexPatchesSelect.get(0)).faces.length / 3 : 0;
        }

        if (keyEvent.getKeyChar() == '\\'){
            stopRendering.set(!stopRendering.get());
        }
        //if (keyEvent.getKeyChar() == 'g'){
        //    convexMeshInitialized = !convexMeshInitialized;
        //    concaveMeshInitialized = !concaveMeshInitialized;
        //    toriMeshInitialized = !toriMeshInitialized;
        //}

        if (keyEvent.getKeyCode() == KeyEvent.VK_F9){
            convexMeshInitialized = !convexMeshInitialized;
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_F10){
            concaveMeshInitialized = !concaveMeshInitialized;
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_F11){
            toriMeshInitialized = !toriMeshInitialized;
        }

        if (Character.isDigit(keyEvent.getKeyChar())){
            strSelectedAtom += keyEvent.getKeyChar();
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_BACK_SPACE){
            strSelectedAtom = "";
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_ENTER){
            try {
                int temp = Integer.parseInt(strSelectedAtom);
                if (temp < convexPatchList.size()){
                    selectedAtom.set(temp);
                } else {
                    strSelectedAtom = "";
                }
            } catch (NumberFormatException e){
                System.out.println("Not a number");
            }
            strSelectedAtom = "";
        }

        if (keyEvent.getKeyChar() == '+'){
            if (selectedAtom.get() + 1 >= convexPatchList.size()){
                selectedAtom.set(0);
            } else {
                selectedAtom.set(selectedAtom.get() + 1);
            }
        }

        if (keyEvent.getKeyChar() == '-'){
            if (selectedAtom.get() - 1 < 0){
                selectedAtom.set(convexPatchList.size() - 1);
            } else {
                selectedAtom.set(selectedAtom.get() - 1);
            }
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_N){
            focusCameraOnTarget();
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_X){
            selectedExclusiveRender = !selectedExclusiveRender;
        }

        if (convexPatchList != null) {
            if (strSelectedAtom.length() > 0) {
                window.setTitle("Selected atom: " + selectedAtom.get() + " / " + convexPatchList.size() + " Atom to select: " + strSelectedAtom + " Press Enter to confirm");
            } else {
                window.setTitle("Selected atom: " + selectedAtom.get() + " / " + convexPatchList.size());
            }
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_ESCAPE){
            captureMouse = !captureMouse;
            rotating = !captureMouse;
            window.confinePointer(captureMouse);
            window.setPointerVisible(!captureMouse);
            window.warpPointer(window.getWidth() / 2, window.getHeight() / 2);
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_W){
            cforward = 1;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_S){
            cforward = -1;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_A){
            cright = -1;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_D){
            cright = 1;
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_K){
            cullFaces = !cullFaces;
            GLRunnable t = (GLAutoDrawable autoDrawable) -> {
                if (cullFaces){
                    gl.glEnable(GL.GL_CULL_FACE);
                } else {
                    gl.glDisable(GL.GL_CULL_FACE);
                }
                return true;
            };
            window.invoke(false, t);
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_F){
            drawMode = (drawMode + 1) % 2; // for now only 'only' modes (face-only, edge-only)
             if (drawMode < 2) {
                GLRunnable t = (GLAutoDrawable drawable) -> {
                    gl.glPolygonMode(GL4.GL_FRONT_AND_BACK, (drawMode == 0) ? GL4.GL_FILL : GL4.GL_LINE);
                    return true;
                };
                window.invoke(false, t);
            }
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_L){
            renderLines = !renderLines;
        }

        if (keyEvent.getKeyChar() == '/'){
            renderCPs = !renderCPs;
        }

        if (keyEvent.getKeyCode() == KeyEvent.VK_SHIFT){
            if (zooming){
                viewPanning = true;
                zooming = false;
            }
            addToSelection = true;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_CONTROL){
            removeSelection = true;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_M){
            mouseSelect = !mouseSelect;
            if (!mouseSelect){
                hoverAtom = -1;
            }
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_SPACE){
            List<Neighbor<double[], SphericalPatch>> neighs = new ArrayList<>();
            SphericalPatch cp = concavePatchList.get(selectedConcaveP.get());
            Surface.probeTree.range(cp.sphere.center.getData(), 2 * Double.longBitsToDouble(Surface.probeRadius.get()), neighs);
            System.out.println("found " + neighs.size() + " neighbors");
            linkNeighbors.clear();
            for (Neighbor<double[], SphericalPatch> n : neighs){
                System.out.println("id: " + n.value.id);
                linkNeighbors.add(n.value.id);
            }
        }
    }
    public static void focusCameraOnTarget(Point p){
        GLRunnable task = new GLRunnable() {
            @Override
            public boolean run(GLAutoDrawable glAutoDrawable) {
                lastCameraTarget[0] = (float)p.x;
                lastCameraTarget[1] = (float)p.y;
                lastCameraTarget[2] = (float)p.z;
                focusCameraOnTarget();
                return true;
            }
        };
       mainWindow.window.invoke(false, task);
    }

    private static void focusCameraOnTarget(){
        camDir = new Quaternion(direction[0], direction[1], direction[2], 0.f);
        camDir.normalize();
        float[] vec = new float[3];
        vec[0] = lastCameraTarget[0] - cameraPos[0];
        vec[1] = lastCameraTarget[1] - cameraPos[1];
        vec[2] = lastCameraTarget[2] - cameraPos[2];
        camTar = new Quaternion(vec[0], vec[1], vec[2], 0.f);
        camTar.normalize();
        slerping = true;
        slerpParam = 0.0f;
    }
    @Override
    public void keyReleased(KeyEvent keyEvent) {
        if (keyEvent.isAutoRepeat()){
            return;
        }
        if (keyEvent.getKeyCode() ==  KeyEvent.VK_W || keyEvent.getKeyCode() == KeyEvent.VK_S){
            cforward = 0;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_A || keyEvent.getKeyCode() == KeyEvent.VK_D){
            cright = 0;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_SHIFT){
            viewPanning = false;
            addToSelection = false;
        }
        if (keyEvent.getKeyCode() == KeyEvent.VK_CONTROL){
            removeSelection = false;
        }
    }
    private int hoverSelectID = -1;

    private int getMeshID(){
        if (hoverAtom > -1){
            if (hoverAtom < convexVerticesCount){
                int localCount = hoverAtom;
                int h = 0;
                for (int i = 0; i < convexPatchList.size(); ++i){
                    SphericalPatch a = convexPatchList.get(i);
                    h += a.vertices.size();
                    if (localCount < h){
                        hoverSelectID = i;
                        break;
                    }
                }
                return CONVEX;
            } else if (hoverAtom < convexVerticesCount + concaveVerticesCount) {
                int localCount = hoverAtom - convexVerticesCount;
                int h = 0;
                for (int i = 0; i < concavePatchList.size(); ++i) {
                    h += concavePatchList.get(i).vertices.size();
                    if (localCount < h) {
                        hoverSelectID = i;
                        break;
                    }
                }
                return CONCAVE;
            } else {
                int localCount = hoverAtom - convexVerticesCount - concaveVerticesCount;
                int h = 0;
                for (int i = 0; i < Surface.rectangles.size(); ++i){
                    if (Surface.rectangles.get(i).faces == null){
                        continue;
                    }
                    h += Surface.rectangles.get(i).faces.length;//.vrts.size() / 2;
                    if (localCount < h){
                        hoverSelectID = i;
                        break;
                    }
                }
                return TORUS;
            }
        }
        return -1;
    }

    @Override
    public void mouseClicked(MouseEvent mouseEvent) {
        if (hoverAtom > -1){
            int meshType = getMeshID();
            if (meshType == CONVEX){
                selectedAtom.set(hoverSelectID);
                if (!addToSelection && !removeSelection){
                    convexPatchesSelect.clear();
                    convexPatchesSelect.add(hoverSelectID);
                } else {
                    if (addToSelection) {
                        convexPatchesSelect.add(hoverSelectID);
                    } else if (removeSelection) {
                        convexPatchesSelect.remove((Object)hoverSelectID);
                    }
                }
            } else if (meshType == CONCAVE){
                selectedConcaveP.set(hoverSelectID);
                if (!addToSelection && !removeSelection){
                    concavePatchesSelect.clear();
                    concavePatchesSelect.add(hoverSelectID);
                } else {
                    if (addToSelection){
                        concavePatchesSelect.add(hoverSelectID);
                    } else if (removeSelection){
                        concavePatchesSelect.remove((Object)hoverSelectID);
                    }
                }
                linkNeighbors.clear();
            } else {
                selectedToriP.set(hoverSelectID);
                if (!addToSelection && !removeSelection){
                    toriPatchesSelect.clear();
                    toriPatchesSelect.add(hoverSelectID);

                } else {
                    if (addToSelection){
                        toriPatchesSelect.add(hoverSelectID);
                    } else if (removeSelection){
                        toriPatchesSelect.remove((Object)hoverSelectID);
                    }
                }
            }
        }
    }

    @Override
    public void mouseEntered(MouseEvent mouseEvent) {

    }

    @Override
    public void mouseExited(MouseEvent mouseEvent) {

    }

    @Override
    public void mousePressed(MouseEvent mouseEvent) {
        if (!mouseDown){
            mouseDown = true;
        }
        if (mouseEvent.getButton() == MouseEvent.BUTTON2){
            if (mouseEvent.isShiftDown()){
                viewPanning = true;
            } else {
                zooming = true;
            }
            zoomSpeed = (mouseEvent.isControlDown()) ? 0.01f : 0.1f;
            window.setPointerVisible(false);
            window.confinePointer(true);
            window.warpPointer(window.getWidth() / 2, window.getHeight() / 2);
            rotating = false;
        }
    }

    @Override
    public void mouseReleased(MouseEvent mouseEvent) {
        mouseDown = false;
        if (zooming) {
            if (!captureMouse) {
                window.setPointerVisible(true);
                window.confinePointer(false);
            }
            zooming = false;
        }
        if (viewPanning){
            if (!captureMouse){
                window.setPointerVisible(true);
                window.confinePointer(false);
            }
            viewPanning = false;
        }
        rotating = true;
    }

    @Override
    public void mouseMoved(MouseEvent mouseEvent) {
        if (captureMouse) {
            mouseAngleX = mouseSpeed * deltaTime * (window.getWidth() / 2.f - mouseEvent.getX());
            mouseAngleY = mouseSpeed * deltaTime * (window.getHeight() / 2.f - mouseEvent.getY());

            if (mouseAngleY > 90.f){
                mouseAngleY = 89.f;
            } else if (mouseAngleY < -90.f){
                mouseAngleY = -89.f;
            }
            axisUp.setFromAngleNormalAxis((float)Math.toRadians(mouseAngleX) * 1.f, up);
            direction = axisUp.rotateVector(direction, 0, direction, 0);
            direction = VectorUtil.normalizeVec3(direction);
            window.warpPointer(window.getWidth() / 2, window.getHeight() / 2);
            if (VectorUtil.dotVec3(direction, new float[]{0.f, 1.f, 0.f}) > Math.cos(Math.toRadians(10)) && mouseAngleY < 0){
                right = VectorUtil.crossVec3(right, new float[]{0.f, 1.f, 0.f}, direction);
                return;
            }
            if (VectorUtil.dotVec3(direction, new float[]{0.f, -1.f, 0.f}) > Math.cos(Math.toRadians(-10)) && mouseAngleY > 0){
                right = VectorUtil.crossVec3(right, new float[]{0.f, 1.f, 0.f}, direction);
                return;
            }
            right = VectorUtil.crossVec3(right, new float[]{0.f, 1.f, 0.f}, direction);

            right = VectorUtil.normalizeVec3(right);
            up = VectorUtil.crossVec3(up, right, direction);
            axisRight.setFromAngleNormalAxis((float)Math.toRadians(mouseAngleY) * 1.f, right);
            direction = axisRight.rotateVector(direction, 0, direction, 0);
            up = VectorUtil.crossVec3(up, right, direction);
            up = VectorUtil.normalizeVec3(up);

        }
        if (mouseSelect){
            mouseLocation.x = mouseEvent.getX();
            mouseLocation.y = mouseEvent.getY();
            moved = true;
        }
    }
    @Override
    public void mouseDragged(MouseEvent mouseEvent) {
        if (zooming){
            if (mouseEvent.isControlDown()){
                zoomSpeed = 0.01f;
            } else {
                zoomSpeed = 0.1f;
            }
            int diffy = (window.getHeight() / 2 - mouseEvent.getY());
            cameraPos[0] += zoomSpeed * direction[0] * diffy;
            cameraPos[1] += zoomSpeed * direction[1] * diffy;
            cameraPos[2] += zoomSpeed * direction[2] * diffy;
            window.warpPointer(window.getWidth() / 2, window.getHeight() / 2);
        }
        if (viewPanning){
            if (mouseEvent.isControlDown()){
                zoomSpeed = 0.01f;
            } else {
                zoomSpeed = 0.1f;
            }
            int diffx = (window.getWidth() / 2 - mouseEvent.getX());
            int diffy = (window.getHeight() / 2 - mouseEvent.getY());
            cameraPos[0] -= zoomSpeed * (-right[0] * diffx + up[0] * diffy);
            cameraPos[1] -= zoomSpeed * (-right[1] * diffx + up[1] * diffy);
            cameraPos[2] -= zoomSpeed * (-right[2] * diffx + up[2] * diffy);
            window.warpPointer(window.getWidth() / 2, window.getHeight() / 2);
        }
    }

    @Override
    public void mouseWheelMoved(MouseEvent mouseEvent) {
        if (mouseEvent.getRotation()[1] < 0.0f){
            scaleFactor -= 0.5f;
        } else if (mouseEvent.getRotation()[1] > 0.0f){
            scaleFactor += 0.5f;
        }
        selectedMeshScaleMat.loadIdentity();
        selectedMeshScaleMat.scale(scaleFactor, scaleFactor, scaleFactor);
    }

    public void freeGLResources(){
        IntBuffer buff = GLBuffers.newDirectIntBuffer(vaos.size());
        convexPatchesSelect.clear();
        concavePatchesSelect.clear();
        convexPatchList = null;
        concavePatchList = null;
        if (vaos.size() > 0) {
            for (Integer i : vaos) {
                buff.put(i);
            }
            buff.rewind();
            gl.glDeleteVertexArrays(vaos.size(), buff);
        }

        if (vbos.size() > 0) {
            buff = GLBuffers.newDirectIntBuffer(vbos.size());
            for (Integer i : vbos) {
                buff.put(i);
            }
            buff.rewind();
            gl.glDeleteBuffers(buff.capacity(), buff);
        }
        if (ebos.size() > 0){
            buff = GLBuffers.newDirectIntBuffer(ebos.size());
            for (Integer i : ebos) {
                buff.put(i);
            }
            buff.rewind();
            gl.glDeleteBuffers(buff.capacity(), buff);
        }
        vaos.clear();
        vbos.clear();
        ebos.clear();
        gl.glDeleteVertexArrays(3, meshVao, 0);
        gl.glDeleteBuffers(3, meshVbo, 0);
        gl.glDeleteBuffers(2, meshEbo, 0);

        gl.glDeleteVertexArrays(2, lineVao, 0);
        gl.glDeleteBuffers(2, lineVbo, 0);
        gl.glDeleteBuffers(2, lineEbo, 0);

        gl.glGenVertexArrays(3, meshVao, 0);
        gl.glGenBuffers(3, meshVbo, 0);
        gl.glGenBuffers(2, meshEbo, 0);

        gl.glGenVertexArrays(2, lineVao, 0);
        gl.glGenBuffers(2, lineVbo, 0);
        gl.glGenBuffers(2, lineEbo, 0);
        resourcesFreed.set(true);
        SurfaceParser.releaseLock();
        selectedExclusiveRender = false;
        if (SesConfig.verbose) {
            System.out.println("GPU resources freed");
        }
    }

    public void close(){
        window.destroy();
    }

    public void setProbeAlpha(float v){
        this.probeAlpha = v;
    }

    public void setMouseSensitivity(float f){
        this.mouseSpeed = f;
    }

    public void stopRendering(boolean v){
        stopRendering.set(v);
    }

    public void requestFreeResources(){
        stopRendering.set(true);
        while (!stoppedRendering.get()){
        }
        GLRunnable task = new GLRunnable() {
            @Override
            public boolean run(GLAutoDrawable glAutoDrawable) {
                freeGLResources();
                return true;
            }
        };
        window.invoke(false, task);
    }

    public void setAmbientStrength(float v){
        this.ambientStrength = v;
    }

    public boolean getResourcesFreed(){
        return resourcesFreed.get();
    }
}

package cz.fi.muni.xmraz3.gui.controllers;

import cz.fi.muni.xmraz3.Surface;
import cz.fi.muni.xmraz3.SesConfig;
import cz.fi.muni.xmraz3.SurfaceParser;
import cz.fi.muni.xmraz3.gui.MainPanel;
import cz.fi.muni.xmraz3.gui.MainWindow;
import javafx.application.Platform;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.concurrent.WorkerStateEvent;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TextField;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.stage.DirectoryChooser;
import javafx.stage.Modality;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class MainPanelController {

    public static DoubleProperty probeRadius = new SimpleDoubleProperty(0.0);
    public static IntegerProperty atomCount = new SimpleIntegerProperty(0);
    @FXML
    private TextField txtConvexPatch;
    @FXML
    private TextField txtConcavePatch;
    @FXML
    private TextField txtRollingPatch;
    @FXML
    private GridPane grdSelectedAtom;
    @FXML
    private Label lblAtomRadius;
    @FXML
    private Label lblBoundaryCount;
    @FXML
    private Label lblAtomName;
    @FXML
    private Button btnOpenFolder;
    @FXML
    private TextField txtFolder;
    @FXML
    private CheckBox chkPinToView;
    @FXML
    private ComboBox<String> cmbResolution;
    @FXML
    private ScrollPane scrOverallInfo;
    @FXML
    GridPane grdOverallInfo;
    @FXML
    Label lblAtomCount;
    @FXML
    private Label lblProbeRadius;
    @FXML
    private CheckBox chkShowProbe;
    @FXML
    private Slider sldProbeAlpha;
    @FXML
    private Slider sldAmbientStrength;
    @FXML
    private ScrollPane scrMainPane;
    @FXML
    private AnchorPane anchScrAnchor;
    @FXML
    private TitledPane tlpExport;
    private static final int tlpExportHeight = 114;
    @FXML
    private TextField txtExportFile;
    @FXML
    private ComboBox<String> cmbExportFormat;
    @FXML
    private Button btnExport;
    @FXML
    private VBox vboSettExp;
    @FXML
    private TitledPane tlpAppSettings;
    private static final int tlpAppSettingsHeight = 224;
    @FXML
    private TitledPane tlpMesh;
    private static final int tlpMeshHeight = 250;
    @FXML
    private Spinner<Double> spinnerEdgeLength;
    @FXML
    private Button btnRemesh;
    private boolean remeshPossible = false;
    @FXML
    private ColorPicker atomColorPick;
    @FXML
    private ColorPicker triangleColorPick;
    @FXML
    private ColorPicker torusColorPick;
    @FXML
    private ColorPicker monoColorPick;
    @FXML
    private Button btnResetColors;
    @FXML
    private CheckBox chkUseMonoColor;
    @FXML
    private Slider sldMouseSensitivity;
    private boolean tlpExclusiveExpand = false;
    public static Stage root;
    public static MainPanelController cont;
    private static final int tbpHeight = 160;
    @FXML
    public void initialize(){
        cont = this;
        Tooltip ttipFolder = new Tooltip("Folder where json files are");
        Tooltip.install(btnOpenFolder, ttipFolder);
        btnOpenFolder.setText("Open");
        btnOpenFolder.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                DirectoryChooser dch = new DirectoryChooser();
                dch.setInitialDirectory(new File(System.getProperty("user.dir")));
                dch.setTitle("Choose folder with json files");
                File selectedFolder = dch.showDialog(root);
                if (selectedFolder == null){
                    return;
                }
                txtFolder.setText(selectedFolder.getAbsolutePath());
                ActionEvent ae = new ActionEvent();
                txtFolder.fireEvent(ae);
            }
        });
        atomColorPick.setValue(Color.rgb(204, 102, 51));
        triangleColorPick.setValue(Color.rgb(31, 143, 0));
        torusColorPick.setValue(Color.rgb(51, 77, 179));
        monoColorPick.setValue(atomColorPick.getValue());
        atomColorPick.valueProperty().addListener(new ChangeListener<Color>() {
            @Override
            public void changed(ObservableValue<? extends Color> observable, Color oldValue, Color newValue) {
                if (MainWindow.mainWindow != null && !chkUseMonoColor.isSelected()) {
                    MainWindow.mainWindow.changeColor(0, (float) newValue.getRed(), (float) newValue.getGreen(), (float) newValue.getBlue());
                }
            }
        });

        triangleColorPick.valueProperty().addListener(new ChangeListener<Color>() {
            @Override
            public void changed(ObservableValue<? extends Color> observable, Color oldValue, Color newValue) {
                if (MainWindow.mainWindow != null && !chkUseMonoColor.isSelected()) {
                    MainWindow.mainWindow.changeColor(1, (float) newValue.getRed(), (float) newValue.getGreen(), (float) newValue.getBlue());
                }
            }
        });

        torusColorPick.valueProperty().addListener(new ChangeListener<Color>() {
            @Override
            public void changed(ObservableValue<? extends Color> observable, Color oldValue, Color newValue) {
                if (MainWindow.mainWindow != null && !chkUseMonoColor.isSelected()) {
                    MainWindow.mainWindow.changeColor(2, (float) newValue.getRed(), (float) newValue.getGreen(), (float) newValue.getBlue());
                }
            }
        });

        monoColorPick.valueProperty().addListener(new ChangeListener<Color>() {
            @Override
            public void changed(ObservableValue<? extends Color> observable, Color oldValue, Color newValue) {
                if (chkUseMonoColor.isSelected()){
                    MainWindow.mainWindow.changeColor(0, (float) monoColorPick.getValue().getRed(), (float) monoColorPick.getValue().getGreen(), (float) monoColorPick.getValue().getBlue());
                    MainWindow.mainWindow.changeColor(1, (float) monoColorPick.getValue().getRed(), (float) monoColorPick.getValue().getGreen(), (float) monoColorPick.getValue().getBlue());
                    MainWindow.mainWindow.changeColor(2, (float) monoColorPick.getValue().getRed(), (float) monoColorPick.getValue().getGreen(), (float) monoColorPick.getValue().getBlue());
                }
            }
        });

        btnResetColors.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                atomColorPick.setValue(Color.rgb(204, 102, 51));
                triangleColorPick.setValue(Color.rgb(31, 143, 0));
                torusColorPick.setValue(Color.rgb(51, 77, 179));
                monoColorPick.setValue(atomColorPick.getValue());
            }
        });
        chkUseMonoColor.selectedProperty().addListener(new ChangeListener<Boolean>() {
            @Override
            public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                if (MainWindow.mainWindow != null) {
                    if (newValue.booleanValue()) {
                        if (MainWindow.mainWindow != null) {
                            MainWindow.mainWindow.changeColor(0, (float) monoColorPick.getValue().getRed(), (float) monoColorPick.getValue().getGreen(), (float) monoColorPick.getValue().getBlue());
                            MainWindow.mainWindow.changeColor(1, (float) monoColorPick.getValue().getRed(), (float) monoColorPick.getValue().getGreen(), (float) monoColorPick.getValue().getBlue());
                            MainWindow.mainWindow.changeColor(2, (float) monoColorPick.getValue().getRed(), (float) monoColorPick.getValue().getGreen(), (float) monoColorPick.getValue().getBlue());
                        }
                    } else {
                        MainWindow.mainWindow.changeColor(0, (float) atomColorPick.getValue().getRed(), (float) atomColorPick.getValue().getGreen(), (float) atomColorPick.getValue().getBlue());
                        MainWindow.mainWindow.changeColor(1, (float) triangleColorPick.getValue().getRed(), (float) triangleColorPick.getValue().getGreen(), (float) triangleColorPick.getValue().getBlue());
                        MainWindow.mainWindow.changeColor(2, (float) torusColorPick.getValue().getRed(), (float) torusColorPick.getValue().getGreen(), (float) torusColorPick.getValue().getBlue());

                    }
                }
            }
        });

        btnExport.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                if (txtExportFile.getLength() == 0){
                    return;
                }
                if (cmbExportFormat.getSelectionModel().getSelectedIndex() == 0){
                    SurfaceParser.exportSTLText(txtExportFile.getText() + ".stl");
                }
                else if (cmbExportFormat.getSelectionModel().getSelectedIndex() == 1){
                    SurfaceParser.exportOBJ(txtExportFile.getText() + ".obj", (char)15);
                }
            }
        });

        btnRemesh.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                btnRemesh.disableProperty().set(true);
                SesConfig.edgeLimit = spinnerEdgeLength.getValue();
                SurfaceParser.remesh();
            }
        });
        txtFolder.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                File f = new File(txtFolder.getText());
                if (!f.exists() || !f.isDirectory()){
                    scrMainPane.requestFocus();
                    txtFolder.setPromptText("No valid directory entered");
                    return;
                }
                if (!checkForValidDirectory(txtFolder.getText())){
                    txtFolder.setText("");
                    txtFolder.setPromptText("No valid directory entered");
                    scrMainPane.requestFocus();
                    return;
                }
                SesConfig.inputFolder = txtFolder.getText();
                startParsingJSON(txtFolder.getText());
            }
        });
        txtConvexPatch.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                try {
                    String[] ids = txtConvexPatch.getText().split(",");
                    List<Integer> ints = new ArrayList<>();
                    for (String id : ids){
                        ints.add(Integer.parseInt(id));
                    }
                    MainPanel.mainView.selectConvexPatchesByIDs(ints);
                } catch (NumberFormatException e){
                    e.printStackTrace();
                }
            }
        });
        txtRollingPatch.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                try {
                    String[] ids = txtRollingPatch.getText().split(",");
                    List<Integer> ints = new ArrayList<>();
                    for (String id : ids){
                        ints.add(Integer.parseInt(id));
                    }
                    MainPanel.mainView.selectToriPatchesByIDs(ints);
                } catch (NumberFormatException e){
                    e.printStackTrace();
                }
            }
        });
        scrMainPane.setHbarPolicy(ScrollPane.ScrollBarPolicy.NEVER);
        scrMainPane.widthProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                anchScrAnchor.setPrefWidth(newValue.doubleValue());
            }
        });
        scrMainPane.heightProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                anchScrAnchor.setPrefHeight(newValue.doubleValue());
            }
        });
        spinnerEdgeLength.valueProperty().addListener(new ChangeListener<Double>() {
            @Override
            public void changed(ObservableValue<? extends Double> observable, Double oldValue, Double newValue) {
                if (!remeshPossible){
                    return;
                }
                if (newValue > 0.01 && newValue < 4){
                    if (Math.abs(SesConfig.edgeLimit - newValue) > 0.0){
                        btnRemesh.disableProperty().set(false);
                    } else {
                        btnRemesh.disableProperty().set(true);
                    }
                }
            }
        });
        SpinnerValueFactory<Double> edgeLenFactory = new SpinnerValueFactory.DoubleSpinnerValueFactory(0.2, 0.9, SesConfig.edgeLimit, 0.05);
        spinnerEdgeLength.setValueFactory(edgeLenFactory);
        SpinnerValueFactory<Double> maxAngleFactory = new SpinnerValueFactory.DoubleSpinnerValueFactory(20, 175, 75, 5);
        txtConcavePatch.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                try {
                    String[] ids = txtConcavePatch.getText().split(",");
                    List<Integer> ints = new ArrayList<>();
                    for (String id : ids){
                        ints.add(Integer.parseInt(id));
                    }
                    MainPanel.mainView.selectConcavePatchesByIDs(ints);
                } catch (NumberFormatException e){
                    e.printStackTrace();
                }
            }
        });
        chkPinToView.setDisable(true);
        chkPinToView.setText("Pin to main view");
        chkPinToView.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                MainPanel.pinnedToView = !MainPanel.pinnedToView;
                Platform.runLater(new Runnable() {
                    @Override
                    public void run() {
                        MainPanel.mainView.isPinned.setValue(MainPanel.pinnedToView);
                    }
                });
            }
        });
        chkShowProbe.setDisable(true);
        chkShowProbe.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                if (MainPanel.mainView != null) {MainPanel.mainView.showProbe(chkShowProbe.isSelected());}
            }
        });
        sldProbeAlpha.setDisable(true);
        sldProbeAlpha.valueProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if (MainPanel.mainView == null){
                    return;
                }
                float v = newValue.intValue() / 100.f;
                MainPanel.mainView.setProbeAlpha(v);
            }
        });
        cmbResolution.getItems().add("800x600");
        cmbResolution.getItems().add("1248x702");
        cmbResolution.getItems().add("1280x720");
        cmbResolution.getItems().add("1920x1080");
        cmbResolution.setDisable(true);
        cmbResolution.getSelectionModel().select(0);
        cmbResolution.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                String[] dims = cmbResolution.getSelectionModel().getSelectedItem().split("x");
                try {
                    if (MainPanel.mainView != null) {
                        MainPanel.mainView.window.setSize(Integer.parseInt(dims[0]), Integer.parseInt(dims[1]));
                    }
                } catch (NumberFormatException e){
                    e.printStackTrace();
                }
            }
        });
        lblProbeRadius.setText("-");
        MainPanel.controlPanel = this;
        probeRadius.addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                lblProbeRadius.setText(Double.toString(newValue.doubleValue()));
                Surface.probeRadius.set(Double.doubleToLongBits(newValue.doubleValue()));
            }
        });
        cmbExportFormat.getItems().add("STL text");
        cmbExportFormat.getItems().add("Wavefront OBJ");
        cmbExportFormat.getSelectionModel().select(0);
        root.setMinWidth(320);
        tlpMesh.setPrefHeight(tlpMeshHeight);
        tlpMesh.setExpanded(true);
        vboSettExp.setPrefHeight(620);
        vboSettExp.setMaxHeight(tlpAppSettingsHeight + tlpMeshHeight + tlpExportHeight + 70 + tbpHeight);
        tlpAppSettings.expandedProperty().addListener(new ChangeListener<Boolean>() {
            @Override
            public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                vboSettExp.setPrefHeight(vboSettExp.getPrefHeight() + ((newValue) ? 1 : -1) * (tlpAppSettingsHeight));
                tlpAppSettings.setPrefHeight((newValue) ? tlpAppSettingsHeight : 0);
                if (tlpExclusiveExpand){
                    tlpExclusiveExpand = false;
                } else {
                    tlpExclusiveExpand = true;
                }
            }
        });
        tlpExport.expandedProperty().addListener(new ChangeListener<Boolean>() {
            @Override
            public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                vboSettExp.setPrefHeight(vboSettExp.getPrefHeight() + ((newValue) ? 1 : -1) * (tlpExportHeight));
                tlpExport.setPrefHeight((newValue) ? tlpExportHeight : 0);
                if (tlpExclusiveExpand){
                    tlpExclusiveExpand = false;
                } else {
                    tlpExclusiveExpand = true;
                }
            }
        });
        tlpMesh.expandedProperty().addListener(new ChangeListener<Boolean>() {
            @Override
            public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                vboSettExp.setPrefHeight(vboSettExp.getPrefHeight() + ((newValue) ? 1 : -1) * tlpMeshHeight);
                tlpMesh.setPrefHeight((newValue) ? tlpMeshHeight: 0);
            }
        });
        sldMouseSensitivity.setDisable(true);
        sldMouseSensitivity.valueProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if (MainPanel.mainView != null){
                    MainPanel.mainView.setMouseSensitivity(newValue.floatValue());
                }
            }
        });
        root.setOnShown(new EventHandler<WindowEvent>() {
            @Override
            public void handle(WindowEvent event) {
                if (SesConfig.inputFolder != null && SesConfig.inputFolder.length() > 0){
                    txtFolder.setText(SesConfig.inputFolder);
                    startParsingJSON(SesConfig.inputFolder);
                }
            }
        });
        sldAmbientStrength.setValue(0.2);
        sldAmbientStrength.setMax(1.0);
        sldAmbientStrength.setMin(0.05);
        sldAmbientStrength.valueProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if (MainWindow.mainWindow != null){
                    MainWindow.mainWindow.setAmbientStrength(newValue.floatValue());
                }
            }
        });
    }

    public void setCheckPinned(boolean v){
        chkPinToView.setSelected(v);
    }

    private void startParsingJSON(String folder){
        btnRemesh.disableProperty().set(true);
        SesConfig.edgeLimit = spinnerEdgeLength.getValue();
        FXMLLoader loader = new FXMLLoader();
        loader.setLocation(MainPanel.class.getResource("/fxml/SESLoadingView.fxml"));
        try {
            String raw = SurfaceParser.loadFile(folder + "/info.json");
            SurfaceParser.parseSesConfig(raw);
            lblAtomCount.setText(Integer.toString(SesConfig.atomCount));
            lblProbeRadius.setText(Double.toString(SesConfig.probeRadius));
            Surface.probeRadius.set(Double.doubleToLongBits(SesConfig.probeRadius));
            AnchorPane anch =  loader.load();
            Scene scene = new Scene(anch);
            Stage stage = new Stage();
            stage.setScene(scene);
            SESLoadingController.stage = stage;
            stage.setTitle("Processing input data");
            SESLoadingController.folder = folder + "/";
            stage.setResizable(false);
            stage.initOwner(root);
            stage.initModality(Modality.APPLICATION_MODAL);
            stage.show();
            SESLoadingController.task.setOnSucceeded(new EventHandler<WorkerStateEvent>() {
                @Override
                public void handle(WorkerStateEvent event) {
                    sldProbeAlpha.setDisable(false);
                    chkShowProbe.setDisable(false);
                    MainPanel.mainView.selectedConcaveP.addListener(new ChangeListener<Number>() {
                        @Override
                        public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                            Platform.runLater(new Runnable() {
                                @Override
                                public void run() {
                                    txtConcavePatch.setText(newValue.toString());
                                }
                            });
                        }
                    });
                    MainPanel.mainView.selectedAtom.addListener(new ChangeListener<Number>() {
                        @Override
                        public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                            Platform.runLater(new Runnable() {
                                @Override
                                public void run() {
                                    txtConvexPatch.setText(newValue.toString());
                                }
                            });
                        }
                    });
                    MainPanel.mainView.selectedToriP.addListener(new ChangeListener<Number>() {
                        @Override
                        public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                            Platform.runLater(new Runnable() {
                                @Override
                                public void run() {
                                    txtRollingPatch.setText(newValue.toString());
                                }
                            });
                        }
                    });
                    sldAmbientStrength.setDisable(false);
                    atomColorPick.setDisable(false);
                    triangleColorPick.setDisable(false);
                    torusColorPick.setDisable(false);
                    monoColorPick.setDisable(false);
                    btnResetColors.setDisable(false);
                    chkUseMonoColor.setDisable(false);
                    MainPanel.mainView.changeColor(0, (float)atomColorPick.getValue().getRed(), (float)atomColorPick.getValue().getGreen(), (float)atomColorPick.getValue().getBlue());
                    MainPanel.mainView.changeColor(1, (float)triangleColorPick.getValue().getRed(), (float)triangleColorPick.getValue().getGreen(), (float)triangleColorPick.getValue().getBlue());
                    MainPanel.mainView.changeColor(2, (float)torusColorPick.getValue().getRed(), (float)torusColorPick.getValue().getGreen(), (float)torusColorPick.getValue().getBlue());
                    chkPinToView.setDisable(false);
                    cmbResolution.setDisable(false);
                    chkPinToView.fire();
                    sldMouseSensitivity.setDisable(false);
                    stage.close();
                }
            });
            SESLoadingController.work();
        } catch (IOException e){
            e.printStackTrace();
        }
    }
    public void updateSelectedAtomInfo(List<String> lst){
        grdSelectedAtom.getChildren().clear();
        for (int i = 0; i < lst.size(); ++i){
            String[] d = lst.get(i).split(",");
            Label key = new Label(d[0] + ":");
            Label value = new Label(d[1]);
            key.setPadding(new Insets(0, 0, 0, 10));
            grdSelectedAtom.addRow(i, key, value);
        }
    }

    private boolean checkForValidDirectory(String dir){
        File f = new File(dir);
        if (!f.isDirectory()){
            return false;
        }
        String[] files = f.list();
        List<String> jsons = new ArrayList<>();
        jsons.add("atoms.dat");
        jsons.add("rectangles.dat");
        jsons.add("triangles.dat");
        jsons.add("info.json");
        int matchCount = 0;
        for (String fileName : files){
            for (String json : jsons){
                if (json.compareTo(fileName) == 0){
                    matchCount++;
                }
            }
        }
        if (matchCount < 4){
            return false;
        }
        return true;
    }

    public static void setBtnRemeshPossible(boolean val){
        cont.remeshPossible = val;
    }

}

package cz.fi.muni.xmraz3.gui;

import cz.fi.muni.xmraz3.gui.controllers.SESLoadingController;
import cz.fi.muni.xmraz3.gui.controllers.MainPanelController;
import javafx.application.Application;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.ObservableList;
import javafx.event.EventHandler;
import javafx.fxml.FXMLLoader;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.layout.AnchorPane;
import javafx.stage.Screen;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;

import java.awt.Rectangle;
import java.awt.GraphicsEnvironment;
import java.io.IOException;

public class MainPanel extends Application {

    public static MainWindow mainView;
    public static boolean pinnedToView = false;
    public static MainPanelController controlPanel;
    @Override
    public void start(Stage primaryStage) {
        MainPanelController.root = primaryStage;
        FXMLLoader loader = new FXMLLoader();
        loader.setLocation(MainPanel.class.getResource("/fxml/MainPanelView.fxml"));
        AnchorPane pane = null;
        try {
            pane = (AnchorPane) loader.load();
            Scene scene = new Scene(pane);
            primaryStage.setScene(scene);

            primaryStage.setOnCloseRequest(new EventHandler<WindowEvent>() {
                @Override
                public void handle(WindowEvent event) {
                    if (mainView != null){
                        mainView.close();
                        mainView = null;
                    }
                    if (SESLoadingController.task != null && SESLoadingController.task.isRunning()){
                        SESLoadingController.task.cancel();
                    }
                    System.exit(0);
                }
            });
            Screen primaryScreen = Screen.getPrimary();
            ObservableList<Screen> screens = Screen.getScreens();
            if (screens.size() > 1){
                Rectangle rect = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getDefaultConfiguration().getBounds();
                for (Screen s : screens){
                    if (s.getBounds().getMinX() == rect.x){
                        primaryScreen = s;
                        break;
                    }
                }
            }
            Rectangle2D primRect = primaryScreen.getVisualBounds();
            primaryStage.setX(primRect.getMinX() + primRect.getWidth() / 2 - scene.getWidth() / 2);
            primaryStage.setY(primRect.getMinY() + primRect.getHeight() / 2 - scene.getHeight() / 2);
            primaryStage.show();
            scene.getWindow().centerOnScreen();
            primaryStage.xProperty().addListener(new ChangeListener<Number>() {
                @Override
                public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                    if (mainView != null && pinnedToView && primaryStage.isFocused()){
                        unpinView();
                    }
                }
            });

            primaryStage.yProperty().addListener(new ChangeListener<Number>() {
                @Override
                public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                    if (mainView != null && pinnedToView && primaryStage.isFocused()){
                        unpinView();
                    }
                }
            });
        } catch (IOException e){
            e.printStackTrace();
        }
        primaryStage.show();
    }

    private void unpinView(){
        pinnedToView = false;
        controlPanel.setCheckPinned(false);
        MainPanel.mainView.isPinned.setValue(false);
    }
    public static void startGUI(String[] args){
        launch(args);
    }
}

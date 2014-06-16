import controlP5.*;
import processing.opengl.*;
import java.io.File;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Properties;
import java.util.StringTokenizer;

ControlP5 cp5;

int natoms = 1;
int nbonds = 1;
int ntimesteps = 1;
int tstep = 0;
boolean play = true;
boolean zoom = false;
int selectedatom = -1;
float atomsize = 0.1;

int framesPerStep = 1;
int framesSinceStep = 0;
    
float[][] x = null;
float[][] y = null;
int[][] bonds = null;
int[] types = null;

float boxwidth = 35.85686;

void setup() {
  size(1100, 800, OPENGL);
  background(30);
  smooth();
  
  readFile();
  
  cp5 = new ControlP5(this);
  cp5.addToggle("play").setPosition(30, 40).setSize(80, 20).setMode(ControlP5.SWITCH).setCaptionLabel("PLAY/PAUSE");
  cp5.addButton("readPositions").setPosition(120, 40).setSize(50, 20).setCaptionLabel("NEXT");
  cp5.addButton("resetPositions").setPosition(180, 40).setSize(50, 20).setCaptionLabel("RESET");
  cp5.addSlider("tstep").setPosition(30, 80).setRange(0, ntimesteps).setSize(200, 20).setCaptionLabel("TIMESTEP");
  cp5.addSlider("framesPerStep").setPosition(30, 120).setRange(1, 20).setSize(200, 20);
}

public void draw() {
  if (play && framesSinceStep++ >= framesPerStep) {
    tstep++;
    framesSinceStep = 1;
  }
  
  if (tstep >= ntimesteps) {
    tstep = 0;
  }
  
  colorMode(RGB, 255);
  background(30);
  
  cp5.controller("tstep").changeValue(tstep);
  
  pushMatrix();
  translate(300, 0);
  scale((width-300)/boxwidth, height/boxwidth);
  
  pushMatrix();
  if (zoom && selectedatom >= 0) {
    float factor = 1.5;
    translate(boxwidth/2.0 - x[tstep][selectedatom], boxwidth/2.0 - y[tstep][selectedatom]);
    scale(factor, factor);
    translate(-boxwidth/(2.0*factor), -boxwidth/(2.0*factor));
  }
  
  noFill();
  stroke(140);
  if (selectedatom >= 0) {
    beginShape();
    for (int t = 0; t < tstep; t++) {
      stroke(((float)t/tstep)*140 + 30);
      vertex(x[t][selectedatom], y[t][selectedatom]);
    }
    endShape();
  }
  
  stroke(60);
  for (int i = 0; i < nbonds; i++) {
    int a = bonds[0][i];
    int b = bonds[1][i];
    if (dist(x[tstep][a], y[tstep][a], x[tstep][b], y[tstep][b]) > 0.5*boxwidth) {
      float xcycle = 0.0;
      float ycycle = 0.0;
      
      if (abs(x[tstep][a] - x[tstep][b]) > 0.5*boxwidth) xcycle = Math.copySign(boxwidth, x[tstep][a] - x[tstep][b]);
      if (abs(y[tstep][a] - y[tstep][b]) > 0.5*boxwidth) ycycle = Math.copySign(boxwidth, y[tstep][a] - y[tstep][b]);
      
      line(x[tstep][a], y[tstep][a], (x[tstep][b] + xcycle), (y[tstep][b] + ycycle));
           
      line((x[tstep][a] - xcycle), (y[tstep][a] - ycycle), x[tstep][b], y[tstep][b]);
    } else {
      line(x[tstep][a], y[tstep][a], x[tstep][b], y[tstep][b]);
    }
  }
  
  colorMode(HSB, 100);
  for (int i = 0; i < natoms; i++) {
    stroke(types[i] * 20, 100, 100);
    fill(types[i] * 20, 100, 100);
    ellipse(x[tstep][i], y[tstep][i], atomsize, atomsize);
  }
  
  popMatrix();
  popMatrix();
}

public void keyPressed() {
  switch (key) {
    case 'h':
      if (cp5.isVisible()) {
        cp5.hide();
      } else {
        cp5.show();
      }
      break;
    case 'n':
      readPositions();
      break;
    case 'r':
      resetPositions();
      break;
    case 's':
      zoom = !zoom;
      break;
    case '-':
      if (atomsize > 0.0) atomsize -= 0.1;
      break;
    case '=':
      atomsize += 0.1;
      break;
  }
}

public void mouseClicked() {
  if (zoom) return;
  
  float mx, my, dx, dy, distSq, minDistSq;
  
  mx = boxwidth*(mouseX - 300)/(width - 300);
  my = boxwidth*mouseY/height;
  minDistSq = Float.POSITIVE_INFINITY;
  
  if (mx < 0 || mx > boxwidth || my < 0 || my > boxwidth) {
    selectedatom = -1;
    return;
  }
  
  for (int i = 0; i < natoms; i++) {
    dx = mx - x[tstep][i];
    dy = my - y[tstep][i];
    distSq = dx*dx + dy*dy;
    
    if (distSq < minDistSq) {
      minDistSq = distSq;
      selectedatom = i;
    }
  }
}

private void readFile() {
  try {
    Properties options = new Properties();
  
    File f = new File(sketchPath, "visualizer.properties");
    options.load(new FileInputStream(f));        
    String directory = options.getProperty("directory");
    File datafile = new File(sketchPath, directory + options.getProperty("datafile"));
    File dumpfile = new File(sketchPath, directory + options.getProperty("dumpfile"));
    File dumpformatfile = new File(sketchPath, directory + options.getProperty("dumpformatfile"));
      
    BufferedReader reader;
    String line;
    
    reader = createReader(datafile);
    reader.readLine();
    reader.readLine();
    line = reader.readLine();
    StringTokenizer st = new StringTokenizer(line.trim());
    natoms = Integer.parseInt(st.nextToken());
    line = reader.readLine();
    st = new StringTokenizer(line.trim());
    nbonds = Integer.parseInt(st.nextToken());
    
    for (int i = 0; i < 10; i++) reader.readLine();
    
    line = reader.readLine();
    st = new StringTokenizer(line.trim());
    st.nextToken();
    boxwidth = Float.parseFloat(st.nextToken());
    
    types = new int[natoms];
    bonds = new int[2][nbonds];
    
    while (!reader.readLine().equals(" Atoms"));
    reader.readLine();
    
    for (int i = 0; i < natoms; i++) {
      line = reader.readLine();
      
      st = new StringTokenizer(line.trim());
      
      if (st.countTokens() != 6) {
        println("Wrong number of tokens in line: " + line);
      } else {
        int id = Integer.parseInt(st.nextToken()) - 1;
        st.nextToken();
        types[id] = Integer.parseInt(st.nextToken());
      }
    }
    
    while (!reader.readLine().equals(" Bonds"));
    reader.readLine();
    
    for (int i = 0; i < nbonds; i++) {
      line = reader.readLine();
      
      st = new StringTokenizer(line.trim());
      
      if (st.countTokens() != 4) {
        println("Wrong number of tokens in line: " + line);
      } else {
        int id = Integer.parseInt(st.nextToken()) - 1;
        st.nextToken();
        bonds[0][id] = Integer.parseInt(st.nextToken()) - 1;
        bonds[1][id] = Integer.parseInt(st.nextToken()) - 1;
      }
    }
    
    reader.close();
    
    // Read parameters from dump format file
    reader = createReader(dumpformatfile);
    reader.readLine(); // Skip text header
    line = reader.readLine();
    reader.close();
    
    st = new StringTokenizer(line); 
    natoms = Integer.parseInt(st.nextToken());
    st.nextToken();
    ntimesteps = Integer.parseInt(st.nextToken());
    
    x = new float[ntimesteps][natoms];
    y = new float[ntimesteps][natoms];
    
    // Read positions from dump file
    FileInputStream dumpstream = new FileInputStream(dumpfile);
    FileChannel dsChannel = dumpstream.getChannel();
    ByteBuffer dsBuffer = ByteBuffer.allocateDirect(8*2*natoms);
    dsBuffer.order(ByteOrder.LITTLE_ENDIAN);
    int bytesRead;
    int t = 0;
    while (t < ntimesteps && (bytesRead = dsChannel.read(dsBuffer)) != -1) {
      dsBuffer.rewind();
      dsBuffer.limit(bytesRead);
      for (int i = 0; i < natoms; i++) {
        x[t][i] = (float) dsBuffer.getDouble();
        y[t][i] = (float) dsBuffer.getDouble();
      }
      t++;
      dsBuffer.clear( );
    }
    
    dumpstream.close();
  } catch (IOException e) {
    e.printStackTrace();
  }
  
  println("Atoms: " + natoms + " Bonds: " + nbonds + " Steps: " + ntimesteps + " Box: " + boxwidth);
}

private void readPositions() {
  tstep++;
}

private void resetPositions() {
  tstep = 0;
}
      

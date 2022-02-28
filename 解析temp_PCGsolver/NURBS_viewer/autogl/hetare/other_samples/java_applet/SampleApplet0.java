import java.awt.*;
import java.awt.event.*;
import java.applet.*;
import java.io.*;

/*
  <applet code="SampleApplet0.class" width=400 height=400>
  </applet>
 */
public class SampleApplet0 extends Applet {

    private void drawView (Graphics g) {
	g.setColor (Color.red);
	g.drawLine (10, 10, 
		    200, 300);

	g.setColor (Color.blue);
	g.drawString (" hello world ",
		      10, 10);
    }



    /* overwride */
    public void init () {
	System.err.println ("init called");

    }

    /* overwride */
    public void paint (Graphics g) {
	System.err.println ("paint called");

	drawView (g);
    }
}


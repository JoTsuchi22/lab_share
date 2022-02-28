import java.awt.*;
import java.awt.event.*;
import java.applet.*;
import java.io.*;

/*
  <applet code="SampleApplet1.class" width=400 height=400>
  </applet>
 */
public class SampleApplet1 extends Applet 
    implements ActionListener {

    private Frame _panelFrame;
    private Button _button0;
    private Button _button1;
    private Checkbox _checkbox0;
    private Checkbox _checkbox1;
    private TextField _textField0;
    private TextField _textField1;
	
    private int _lineStartX = 100;
    private int _lineStartY = 50;
    private int _lineEndX = 200;
    private int _lineEndY = 250;

    private void makePanel () {
	_panelFrame = new Frame ();
	_panelFrame.setTitle ("sample1");
	_panelFrame.setLocation (100, 100);
	_panelFrame.setSize (300, 800);
	_panelFrame.setResizable (false);
	_panelFrame.setVisible (true);

	_panelFrame.setLayout (new GridLayout (0, 1));

	_panelFrame.add (new Label ("label"));

	_button0 = new Button ("button 0");
	_button0.addActionListener (this);
	_panelFrame.add (_button0);

	_button1 = new Button ("button 1");
	_button1.addActionListener (this);
	_panelFrame.add (_button1);

	_checkbox0 = new Checkbox ("checkbox 0");
	_panelFrame.add (_checkbox0);

	_checkbox1 = new Checkbox ("checkbox 1");
	_panelFrame.add (_checkbox1);

	_panelFrame.add (new Label ("text field 0"));
	_textField0 = new TextField (20);
	_panelFrame.add (_textField0);

	_panelFrame.add (new Label ("text field 1"));
	_textField1 = new TextField (20);
	_panelFrame.add (_textField1);

	_panelFrame.setResizable (true);
	_panelFrame.pack ();
    }

    private void performButtonAction0 () {
	/* do nothing */
    }

    private void performButtonAction1 () {
	_lineStartX += 10;
	_lineEndX += 10;
    }

    private void performDrag 
	(int startX, int startY, int endX, int endY) {
	_lineStartX = startX;
	_lineStartY = startY;
	_lineEndX = endX;
	_lineEndY = endY;
    }

    private void drawView (Graphics g) {
	if (_checkbox0.getState ()) {
	    g.setColor (Color.red);
	} else {
	    g.setColor (Color.black);
	}
	g.drawLine (_lineStartX, _lineStartY, 
		    _lineEndX, _lineEndY);


	if (_checkbox1.getState ()) {
	    g.setColor (Color.red);
	} else {
	    g.setColor (Color.black);
	}
	g.drawString (_textField0.getText (),
		      _lineStartX, _lineStartY);
	g.drawString (_textField1.getText (),
		      _lineEndX, _lineEndY);
    }



    /* overwride */
    public void init () {
	System.err.println ("init called");

	makePanel ();
    }

    /* overwride */
    public void paint (Graphics g) {
	System.err.println ("paint called");

	drawView (g);
    }
    
    /* overwride */
    public void actionPerformed (ActionEvent e) {
	System.err.println ("actionPerformed called");

	Button source = (Button)e.getSource ();
	if (source == _button0) {
	    performButtonAction0 ();

	    repaint ();
	} else if (source == _button1) {
	    performButtonAction1 ();

	    repaint ();
	}
    }

    private int _startX;
    private int _startY;

    /* overwride */
    public boolean mouseDown (Event e, int x, int y) {
	System.err.println ("mouseDown called");

	_startX = x;
	_startY = y;
	return true;
    }

    /* overwride */
    public boolean mouseDrag (Event e, int x, int y) {
	System.err.println ("mouseDrag called");

	return true;
    }

    /* overwride */
    public boolean mouseUp (Event e, int x, int y) {
	System.err.println ("mouseUp called");

	performDrag (_startX, _startY, x, y);

	repaint ();
	return true;
    }
}


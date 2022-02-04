import java.awt.*;
import java.awt.event.*;
import java.applet.*;
import java.io.*;
import java.net.*;
import java.util.*;

/*
  <applet code="AutoGLApplet.class" width=500 height=500>
  </applet>
 */
public class AutoGLApplet extends Applet 
    implements ActionListener {

    /* basic */

    private void error (String message) {
	System.err.println ("ERROR");
	System.err.println (message);
    }



    /* reader */

    private static final int MAX_WORD_CHARS = 256;

    private InputStreamReader _reader;
    private int _currentChar;
    private char _wordChars[] = new char[MAX_WORD_CHARS];

    private void initializeReader (InputStreamReader reader) {
	_reader = reader;
	try {
	    _currentChar = _reader.read ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void terminateReader () {
	try {
	    _reader.close ();
	    _reader = null;
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private boolean isEOF () {
	return _currentChar == -1;
    }

    private boolean isWhiteSpace () {
	return _currentChar == ' ' 
	    || _currentChar == '\n'
	    || _currentChar == '\t';
    }

    private String readWord () {
	try {
	    while (isWhiteSpace ()) {
		_currentChar = _reader.read ();
	    }
	    if (isEOF ()) return null;
	    
	    int count = 0;
	    if (_currentChar == '"') {
		_currentChar = _reader.read ();
		for (;;) {
		    if (isEOF ()) {
			error ("must be \"");
			break;
		    }
		    if (_currentChar == '"') {
			_currentChar = _reader.read ();
			break;
		    }
		    
		    if (MAX_WORD_CHARS <= count) {
			error ("word is too long");
			break;
		    }
		    _wordChars[count++] = (char)_currentChar;
		    _currentChar = _reader.read ();
		}
	    } else {
		for (;;) {
		    if (isWhiteSpace ()) break;
		    if (isEOF ()) break;
		    
		    if (MAX_WORD_CHARS <= count) {
			error ("word is too long");
			break;
		    }
		    _wordChars[count++] = (char)_currentChar;
		    _currentChar = _reader.read ();
		}
	    }
	    return new String (_wordChars, 0, count);
	} catch (IOException ex) {
	    error (ex.toString ());
	}
	return null;
    }
    


    /* view & panel */

    private Frame _panelFrameComponent;

    private void makePanel () {
	_panelFrameComponent = new Frame ();
	_panelFrameComponent.setTitle ("AutoGL:Panel");
	_panelFrameComponent.setLocation (100, 100);
	_panelFrameComponent.setSize (300, 800);
	_panelFrameComponent.setResizable (false);
	_panelFrameComponent.setVisible (true);

	_panelFrameComponent.setLayout (new GridLayout (0, 1));

	readModelLayout ();
	readModelValue ();

	_panelFrameComponent.setResizable (true);
	_panelFrameComponent.pack ();
    }

    private void makeView () {
	postDrawCommand ();
	readImage ();
    }



    /* widget */

    private static final int MAX_WIDGETS = 100;
    private static final int MAX_LABEL_COMPONENTS = 100;
    private static final int MAX_COMPONENTS = 50;

    private static final int PROPERTY_GROUP = 0;
    private static final int PROPERTY_COMMENT = 1;
    private static final int PROPERTY_CALLBACK = 2;
    private static final int PROPERTY_BOOLEAN = 3;
    private static final int PROPERTY_STRING = 4;
    private static final int PROPERTY_INTEGER = 5;
    private static final int PROPERTY_REAL = 6;

    private static final int WIDGET_LABEL = 0;
    private static final int WIDGET_BUTTON = 1;
    private static final int WIDGET_TOGGLE = 2;
    private static final int WIDGET_FIELD = 3;
    private static final int WIDGET_CHOICE = 4;

    private int _nWidgets;
    private String _widgetNames[];
    private int _widgetTypes[];
    private int _widgetComponentIndexs[];
    private int _widgetPropertyTypes[];

    private int _nButtonComponents;
    private Button _buttonComponents[];
    private int _buttonWidgetIds[];

    private int _nCheckboxComponents;
    private Checkbox _checkboxComponents[];

    private int _nTextFieldComponents;
    private TextField _textFieldComponents[];

    private int _nChoiceComponents;
    private Choice _choiceComponents[];

    private void clearWidget () {
	_nWidgets = 0;
	_widgetNames = new String[MAX_WIDGETS];
	_widgetTypes = new int[MAX_WIDGETS];
	_widgetComponentIndexs = new int[MAX_WIDGETS];
	_widgetPropertyTypes = new int[MAX_WIDGETS];
	
	_nButtonComponents = 0;
	_buttonComponents = new Button[MAX_COMPONENTS];
	_buttonWidgetIds = new int[MAX_COMPONENTS];
	
	_nCheckboxComponents = 0;
	_checkboxComponents = new Checkbox[MAX_COMPONENTS];
	
	_nTextFieldComponents = 0;
	_textFieldComponents = new TextField[MAX_COMPONENTS];

	_nChoiceComponents = 0;
	_choiceComponents = new Choice[MAX_COMPONENTS];
    }

    private int widgetIdOfName (String name) {
	for (int iWidget = 0; iWidget < _nWidgets; iWidget++) {
	    if (_widgetNames[iWidget].equals (name)) {
		return iWidget;
	    }
	}
	return -1;
    }



    /* label widget */

    private void makeLabelWidget 
	(String name, String label, int propertyType) {
	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_LABEL;
	_widgetComponentIndexs[_nWidgets] = -1;
	_widgetPropertyTypes[_nWidgets] = propertyType;
	_nWidgets++;

	_panelFrameComponent.add (new Label (label));
    }



    /* button widget */

    private void makeButtonWidget (String name, String label) {
	int widgetId = _nWidgets;
	int componentIndex = _nButtonComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_BUTTON;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = PROPERTY_CALLBACK;
	_nWidgets++;

	String buttonLabel;
	if (0 < label.length ()) {
	    buttonLabel = label;
	} else {
	    buttonLabel = name;
	}
	Button buttonComponent = new Button (buttonLabel);
	buttonComponent.addActionListener (this);
	_panelFrameComponent.add (buttonComponent);

	_buttonComponents[_nButtonComponents] = buttonComponent;
	_buttonWidgetIds[_nButtonComponents] = widgetId;
	_nButtonComponents++;
    }



    /* toggle widget */

    private void makeToggleWidget (String name, String label) {
	int componentIndex = _nCheckboxComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_TOGGLE;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = PROPERTY_BOOLEAN;
	_nWidgets++;
	
	String checkboxLabel;
	if (0 < label.length ()) {
	    checkboxLabel = label;
	} else {
	    checkboxLabel = name;
	}
	Checkbox checkboxComponent = new Checkbox (checkboxLabel);
	_panelFrameComponent.add (checkboxComponent);

	_checkboxComponents[_nCheckboxComponents] = checkboxComponent;
	_nCheckboxComponents++;
    }

    private Checkbox componentOfToggleWidget (int widgetId) {
	if (_widgetTypes[widgetId] != WIDGET_TOGGLE) {
	    error ("type must be toggle");
	}
	int componentIndex = _widgetComponentIndexs[widgetId];
	return _checkboxComponents[componentIndex];
    }



    /* field widget */

    private static final int MAX_FIELD_CHARS = 16;

    private void makeFieldWidget 
	(String name, String label, int propertyType) {
	int componentIndex = _nTextFieldComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_FIELD;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = propertyType;
	_nWidgets++;
	
	Panel panel = new Panel ();
	panel.setLayout (new FlowLayout ());
	
	String fieldLabel;
	if (0 < label.length ()) {
	    fieldLabel = label;
	} else {
	    fieldLabel = name;
	}
	Label labelComponent = new Label (fieldLabel);
	panel.add (labelComponent);
	
	TextField textFieldComponent = new TextField (MAX_FIELD_CHARS);
	panel.add (textFieldComponent);

	_textFieldComponents[_nTextFieldComponents] = textFieldComponent;
	_nTextFieldComponents++;

	_panelFrameComponent.add (panel);
    }

    private TextField componentOfFieldWidget (int widgetId) {
	if (_widgetTypes[widgetId] != WIDGET_FIELD) {
	    error ("type must be field");
	}
	int componentIndex = _widgetComponentIndexs[widgetId];
	return _textFieldComponents[componentIndex];
    }



    /* choice widget */

    private void makeChoiceWidget (String name, String label) {
	int componentIndex = _nChoiceComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_CHOICE;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = PROPERTY_INTEGER;
	_nWidgets++;
	
	Panel panel = new Panel ();
	panel.setLayout (new FlowLayout ());
	
	String choiceLabel;
	if (0 < label.length ()) {
	    choiceLabel = label;
	} else {
	    choiceLabel = name;
	}
	Label labelComponent = new Label (choiceLabel);
	panel.add (labelComponent);

	Choice choiceComponent = new Choice ();
	panel.add (choiceComponent);

	_choiceComponents[_nChoiceComponents] 
	    = choiceComponent;
	_nChoiceComponents++;

	_panelFrameComponent.add (panel);
    }

    private void addItemIntoChoiceWidget (String itemLabel) {
	Choice currentChoiceComponent 
	    = _choiceComponents[_nChoiceComponents - 1];
	currentChoiceComponent.add (itemLabel);
    }

    private Choice componentOfChoiceWidget (int widgetId) {
	if (_widgetTypes[widgetId] != WIDGET_CHOICE) {
	    error ("type must be choice");
	}
	int componentIndex = _widgetComponentIndexs[widgetId];
	return _choiceComponents[componentIndex];
    }



    /* application model I/O */

    private Image _viewImage;

    private String VIEW_IMAGE = "autogl_view_image";
    private String VIEW_IMAGE_COUNT = "autogl_view_image_count_file";
    private String INFO_COMMAND_OUTPUT = "autogl_info_command_output_file";
    private String VALUE_COMMAND_OUTPUT = "autogl_value_command_output_file";
    private String POST_COMMAND_CGI = "autogl_post_command.cgi";

    private static final int VIEW_IMAGE_WIDTH = 500;
    private static final int VIEW_IMAGE_HEIGHT = 500;

    private void readImage () {
	int count = -1;
	try {
	    URL url = new URL (getDocumentBase (), VIEW_IMAGE_COUNT);
	    InputStreamReader r 
		= new InputStreamReader (url.openStream ());
	    BufferedReader br = new BufferedReader (r);
	    String value = br.readLine ();
	    count = Integer.parseInt (value);
	    br.close ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}

	// String fileName = VIEW_IMAGE + count + ".gif";
	String fileName = VIEW_IMAGE + count + ".jpeg";

	_viewImage = getImage (getDocumentBase (), fileName);
    }

    private void parseInfoCommandOutputFile () {
	for (;;) {
	    String keyword = readWord ();
	    if (keyword == null) break;
	
	    String nameKeyword = readWord ();
	    if (!nameKeyword.equals ("name")) {
		error ("must be name");
		return;
	    }
	    String name = readWord ();
	    
	    String labelKeyword = readWord ();
	    if (!labelKeyword.equals ("label")) {
		error ("must be label");
		return;
	    }
	    String label = readWord ();
	    
	    if (keyword.equals ("group")) {
		makeLabelWidget (name, label, PROPERTY_GROUP);
		
	    } else if (keyword.equals ("comment")) {
		makeLabelWidget (name, label, PROPERTY_COMMENT);
		
	    } else if (keyword.equals ("callback")) {
		makeButtonWidget (name, label);
		
	    } else if (keyword.equals ("boolean")) {
		makeToggleWidget (name, label);
		
	    } else if (keyword.equals ("string")) {
		makeFieldWidget (name, label, PROPERTY_STRING);
		
	    } else if (keyword.equals ("integer")) {
		String nItemsKeyword = readWord ();
		if (!nItemsKeyword.equals ("nItems")) {
		    error ("must be nItems");
		    return;
		}
		String nItemsStr = readWord ();
		int nItems = Integer.parseInt (nItemsStr);
		if (nItems == 0) {
		    makeFieldWidget (name, label, PROPERTY_INTEGER);
		} else {
		    makeChoiceWidget (name, label);
		}
		for (int iItem = 0; iItem < nItems; iItem++) {
		    String itemLabelKeyword = readWord ();
		    if (!itemLabelKeyword.equals ("itemLabel")) {
			error ("must be itemLabel");
		    }
		    String itemLabel = readWord ();
		    addItemIntoChoiceWidget (itemLabel);
		}
		
	    } else if (keyword.equals ("real")) {
		makeFieldWidget (name, label, PROPERTY_REAL);
		
	    } else {
		error ("unsupported widget type");
		return;
	    }
	}
    }

    private void readModelLayout () {
	try {
	    URL url = new URL (getDocumentBase (), INFO_COMMAND_OUTPUT);
	    initializeReader 
		(new InputStreamReader (url.openStream ()));
	    parseInfoCommandOutputFile ();
	    terminateReader ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}

	/*
	makeLabelWidget ("LABEL", "label", PROPERTY_COMMENT);
	makeButtonWidget ("BUTTON", "button");
	makeToggleWidget ("TOGGLE", "toggle");
	makeFieldWidget ("FIELD", "field", PROPERTY_STRING);
	makeChoiceWidget ("CHOICE", "choice");
	addItemIntoChoiceWidget ("choice0");
	addItemIntoChoiceWidget ("choice1");
	addItemIntoChoiceWidget ("choice2");
	makeFieldWidget ("REAL_FIELD", "real", PROPERTY_REAL);
	*/
    }

    private void parseValueCommandOutputFile () {
	for (;;) {
	    String keyword = readWord ();
	    if (keyword == null) break;
	
	    String nameKeyword = readWord ();
	    if (!nameKeyword.equals ("name")) {
		error ("must be name");
		return;
	    }
	    String name = readWord ();
	    
	    if (keyword.equals ("group")
		|| keyword.equals ("comment")
		|| keyword.equals ("callback")) continue;
		
	    String valueKeyword = readWord ();
	    if (!valueKeyword.equals ("value")) {
		error ("must be value");
		return;
	    }
	    String value = readWord ();
	    
	    if (keyword.equals ("boolean")) {
		int widgetId = widgetIdOfName (name);
		Checkbox component = componentOfToggleWidget (widgetId);
		if (value.equals ("true")) {
		    component.setState (true);
		} else if (value.equals ("false")) {
		    component.setState (false);
		} else {
		    error ("value must be true or false");
		    return;
		}
		
	    } else if (keyword.equals ("string")) {
		int widgetId = widgetIdOfName (name);
		TextField component = componentOfFieldWidget (widgetId);
		component.setText (value);
		
	    } else if (keyword.equals ("integer")) {
		int widgetId = widgetIdOfName (name);
		if (_widgetTypes[widgetId] == WIDGET_FIELD) {
		    TextField component = componentOfFieldWidget (widgetId);
		    component.setText (value);

		} else if (_widgetTypes[widgetId] == WIDGET_CHOICE) {
		    int integerValue = Integer.parseInt (value);
		    String itemLabel = readWord ();

		    Choice component = componentOfChoiceWidget (widgetId);
		    component.select (integerValue);

		} else {
		    error ("integer must be field or choice");
		    return;
		}
		
	    } else if (keyword.equals ("real")) {
		int widgetId = widgetIdOfName (name);
		TextField component = componentOfFieldWidget (widgetId);
		component.setText (value);
		
	    } else {
		error ("unsupported widget type");
		return;
	    }
	}
    }

    private void readModelValue () {
	try {
	    URL url = new URL (getDocumentBase (), VALUE_COMMAND_OUTPUT);
	    initializeReader 
		(new InputStreamReader (url.openStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void writeResizeCommand (OutputStreamWriter writer) {
	try {
	    writer.write ("resize");
	    writer.write (" ");
	    writer.write (String.valueOf (VIEW_IMAGE_WIDTH));
	    writer.write (" ");
	    writer.write (String.valueOf (VIEW_IMAGE_HEIGHT));
	    writer.write ("\n");
	} catch (Exception ex) {
	    error (ex.toString ());
	}
    }

    private void writeSetCommand 
	(OutputStreamWriter writer,
	 String name, String value) {
	try {
	    writer.write ("set");
	    writer.write (" ");
	    writer.write (name);
	    writer.write (" ");
	    writer.write (value);
	    writer.write ("\n");
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void writeSetCommandForEachWidget 
	(OutputStreamWriter writer) {
	for (int iWidget = 0; iWidget < _nWidgets; iWidget++) {
	    String name = _widgetNames[iWidget];
	    switch (_widgetTypes[iWidget]) {
	    case WIDGET_LABEL:
	    case WIDGET_BUTTON:
		break;
	    case WIDGET_TOGGLE:
		{
		    Checkbox component = 
			componentOfToggleWidget (iWidget);
		    if (component.getState ()) {
			writeSetCommand (writer, name, "true");
		    } else {
			writeSetCommand (writer, name, "false");
		    }
		}
		break;
	    case WIDGET_FIELD:
		{
		    TextField component = 
			componentOfFieldWidget (iWidget);
		    switch (_widgetPropertyTypes[iWidget]) {
		    case PROPERTY_STRING:
			writeSetCommand 
			    (writer, name, 
			     "\"" + component.getText () + "\"");
			break;
		    case PROPERTY_INTEGER:
		    case PROPERTY_REAL:
			writeSetCommand (writer, name, component.getText ());
			break;
		    default:
			error ("illegal property type of field");
			break;
		    }
		}
		break;
	    case WIDGET_CHOICE:
		{
		    Choice component = 
			componentOfChoiceWidget (iWidget);
		    String value = component.getSelectedItem ();
		    writeSetCommand (writer, name, value);
		}
		break;
	    default:
		error ("illegal widget type");
		break;
	    }
	}
    }
    
    private void writeCallCommand 
	(OutputStreamWriter writer, String callbackName) {
	try {
	    writer.write ("call");
	    writer.write (" ");
	    writer.write (callbackName);
	    writer.write ("\n");
	} catch (Exception ex) {
	    error (ex.toString ());
	}
    }
    
    private void writeDrawCommand (OutputStreamWriter writer) {
	try {
	    writer.write ("draw");
	    writer.write ("\n");
	} catch (Exception ex) {
	    error (ex.toString ());
	}
    }

    private void writeDragCommand 
	(OutputStreamWriter writer,
	 int startX, int startY, int endX, int endY) {
	try {
	    writer.write ("drag");
	    writer.write (" ");
	    writer.write (String.valueOf (startX));
	    writer.write (" ");
	    writer.write (String.valueOf (startY));
	    writer.write (" ");
	    writer.write (String.valueOf (endX));
	    writer.write (" ");
	    writer.write (String.valueOf (endY));
	    writer.write ("\n");
	} catch (Exception ex) {
	    error (ex.toString ());
	}
    }

    private void writeValueCommand (OutputStreamWriter writer) {
	try {
	    writer.write ("value");
	    writer.write ("\n");
	} catch (Exception ex) {
	    error (ex.toString ());
	}
    }

    private void postCallCommand (String callbackName) {
	try {
	    URL url = new URL (getDocumentBase (), POST_COMMAND_CGI);
	    URLConnection connection = url.openConnection ();
	    connection.setDoInput (true);
	    connection.setDoOutput (true);
	    connection.setUseCaches (false);

	    OutputStreamWriter writer 
		= new OutputStreamWriter (connection.getOutputStream ());
	    writeResizeCommand (writer);
	    writeSetCommandForEachWidget (writer);
	    writeCallCommand (writer, callbackName);
	    writeValueCommand (writer);
	    writer.close ();

	    initializeReader 
		(new InputStreamReader (connection.getInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void postDrawCommand () {
	try {
	    URL url = new URL (getDocumentBase (), POST_COMMAND_CGI);
	    URLConnection connection = url.openConnection ();
	    connection.setDoInput (true);
	    connection.setDoOutput (true);
	    connection.setUseCaches (false);

	    OutputStreamWriter writer 
		= new OutputStreamWriter (connection.getOutputStream ());
	    writeResizeCommand (writer);
	    writeSetCommandForEachWidget (writer);
	    writeDrawCommand (writer);
	    writeValueCommand (writer);
	    writer.close ();

	    initializeReader 
		(new InputStreamReader (connection.getInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void postDragCommand 
	(int startX, int startY, int endX, int endY) {
	try {
	    URL url = new URL (getDocumentBase (), POST_COMMAND_CGI);
	    URLConnection connection = url.openConnection ();
	    connection.setDoInput (true);
	    connection.setDoOutput (true);
	    connection.setUseCaches (false);

	    OutputStreamWriter writer 
		= new OutputStreamWriter (connection.getOutputStream ());
	    writeResizeCommand (writer);
	    writeSetCommandForEachWidget (writer);
	    writeDragCommand (writer, startX, startY, endX, endY);
	    writeDrawCommand (writer);
	    writeValueCommand (writer);
	    writer.close ();

	    initializeReader 
		(new InputStreamReader (connection.getInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    /* overwride */
    public void init () {
	clearWidget ();
	makePanel ();
	makeView ();
    }

    /* overwride */
    public void paint (Graphics g) {
	g.drawImage (_viewImage, 0, 0, this);
    }
    
    /* overwride */
    public void actionPerformed (ActionEvent e) {
	Button source = (Button)e.getSource ();
	for (int iComponent = 0; 
	     iComponent < _nButtonComponents; 
	     iComponent++) {
	    if (source == _buttonComponents[iComponent]) {
		int widgetId = _buttonWidgetIds[iComponent];
		String name = _widgetNames[widgetId];
		postCallCommand (name);
		readImage ();
		repaint ();
	    }
	}
    }

    private int _startX;
    private int _startY;

    /* overwride */
    public boolean mouseDown (Event e, int x, int y) {
	_startX = x;
	_startY = y;
	return true;
    }

    /* overwride */
    public boolean mouseDrag (Event e, int x, int y) {
	return true;
    }

    /* overwride */
    public boolean mouseUp (Event e, int x, int y) {
	postDragCommand (_startX, _startY, x, y);
	readImage ();
	repaint ();
	return true;
    }
}


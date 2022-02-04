import java.io.*;
import javax.microedition.io.*;
import com.nttdocomo.io.*;
import com.nttdocomo.ui.*;

class AutoGLViewCanvas extends Canvas {

    private AutoGLApplication _application;
    
    private static final int DRAG_CENTER_X = 120;
    private static final int DRAG_CENTER_Y = 120;
    private static final int DRAG_SIZE_X = 60;
    private static final int DRAG_SIZE_Y = 60;

    public void initialize (AutoGLApplication application) {
	_application = application;
	setSoftLabel (Frame.SOFT_KEY_2, "Panel");
    }

    /* overwride */
    public void paint (Graphics g) {
	g.lock ();
	g.clearRect (0, 0, getWidth (), getHeight ());
	g.drawImage (_application.viewImage (), 0, 0);
	g.unlock (true);
    }
    
    /* overwride */
    public void processEvent (int type, int param) {
	if (type == Display.KEY_RELEASED_EVENT) {
	    if (param == Display.KEY_SOFT2) {
		_application.switchToPanel ();
	    } else if (param == Display.KEY_SELECT) {
		_application.drawView ();
		repaint ();
	    } else if (param == Display.KEY_1) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X - DRAG_SIZE_X,
		     DRAG_CENTER_Y - DRAG_SIZE_Y);
		repaint ();
	    } else if (param == Display.KEY_2) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X,
		     DRAG_CENTER_Y - DRAG_SIZE_Y);
		repaint ();
	    } else if (param == Display.KEY_3) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X + DRAG_SIZE_X,
		     DRAG_CENTER_Y - DRAG_SIZE_Y);
		repaint ();
	    } else if (param == Display.KEY_4) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X - DRAG_SIZE_X,
		     DRAG_CENTER_Y);
		repaint ();
	    } else if (param == Display.KEY_5) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X,
		     DRAG_CENTER_Y);
		repaint ();
	    } else if (param == Display.KEY_6) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X + DRAG_SIZE_X,
		     DRAG_CENTER_Y);
		repaint ();
	    } else if (param == Display.KEY_7) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X - DRAG_SIZE_X,
		     DRAG_CENTER_X + DRAG_SIZE_Y);
		repaint ();
	    } else if (param == Display.KEY_8) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X,
		     DRAG_CENTER_Y + DRAG_SIZE_Y);
		repaint ();
	    } else if (param == Display.KEY_9) {
		_application.dragView 
		    (DRAG_CENTER_X, 
		     DRAG_CENTER_Y, 
		     DRAG_CENTER_X + DRAG_SIZE_X,
		     DRAG_CENTER_Y + DRAG_SIZE_Y);
		repaint ();
	    }
	}
    }
} 

public class AutoGLApplication extends IApplication 
    implements ComponentListener, SoftKeyListener {

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

    private Panel _panelPanelComponent;
    private HTMLLayout _panelHTMLLayout;
    private AutoGLViewCanvas _viewCanvasComponent;

    private void makePanel () {
	_panelPanelComponent = new Panel ();
	_panelHTMLLayout = new HTMLLayout ();
	_panelPanelComponent.setLayoutManager (_panelHTMLLayout);
	_panelPanelComponent.setComponentListener (this);
	_panelPanelComponent.setSoftKeyListener (this);
	_panelPanelComponent.setSoftLabel (Frame.SOFT_KEY_2, "View");

	_panelHTMLLayout.begin (HTMLLayout.LEFT);
	readModelLayout ();
	readModelValue ();
	_panelHTMLLayout.end ();
    }

    private void makeView () {
	_viewCanvasComponent = new AutoGLViewCanvas ();
	_viewCanvasComponent.initialize (this);

	postDrawCommand ();
	readImage ();
    }

    private void switchToView () {
	Display.setCurrent (_viewCanvasComponent);
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

    private int _nCheckBoxComponents;
    private ListBox _checkBoxComponents[];

    private int _nTextBoxComponents;
    private TextBox _textBoxComponents[];

    private int _nChoiceComponents;
    private ListBox _choiceComponents[];

    private void clearWidget () {
	_nWidgets = 0;
	_widgetNames = new String[MAX_WIDGETS];
	_widgetTypes = new int[MAX_WIDGETS];
	_widgetComponentIndexs = new int[MAX_WIDGETS];
	_widgetPropertyTypes = new int[MAX_WIDGETS];
	
	_nButtonComponents = 0;
	_buttonComponents = new Button[MAX_COMPONENTS];
	_buttonWidgetIds = new int[MAX_COMPONENTS];
	
	_nCheckBoxComponents = 0;
	_checkBoxComponents = new ListBox[MAX_COMPONENTS];
	
	_nTextBoxComponents = 0;
	_textBoxComponents = new TextBox[MAX_COMPONENTS];

	_nChoiceComponents = 0;
	_choiceComponents = new ListBox[MAX_COMPONENTS];
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

	_panelPanelComponent.add (new Label (label));

	_panelHTMLLayout.br ();
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

	Button buttonComponent;
	if (0 < label.length ()) {
	    buttonComponent = new Button (label);
	} else {
	    buttonComponent = new Button (name);
	}
	_panelPanelComponent.add (buttonComponent);
	_buttonComponents[_nButtonComponents] = buttonComponent;
	_buttonWidgetIds[_nButtonComponents] = widgetId;
	_nButtonComponents++;

	_panelHTMLLayout.br ();
    }



    /* toggle widget */

    private void makeToggleWidget (String name, String label) {
	int componentIndex = _nCheckBoxComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_TOGGLE;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = PROPERTY_BOOLEAN;
	_nWidgets++;
	
	ListBox checkBoxComponent = new ListBox (ListBox.CHECK_BOX);
	if (0 < label.length ()) {
	    checkBoxComponent.append (label);
	} else {
	    checkBoxComponent.append (name);
	}
	_panelPanelComponent.add (checkBoxComponent);
	_checkBoxComponents[_nCheckBoxComponents] = checkBoxComponent;
	_nCheckBoxComponents++;

	_panelHTMLLayout.br ();
    }

    private ListBox componentOfToggleWidget (int widgetId) {
	if (_widgetTypes[widgetId] != WIDGET_TOGGLE) {
	    error ("type must be toggle");
	}
	int componentIndex = _widgetComponentIndexs[widgetId];
	return _checkBoxComponents[componentIndex];
    }



    /* field widget */

    private static final int MAX_FIELD_CHARS = 16;

    private void makeFieldWidget 
	(String name, String label, int propertyType) {
	int componentIndex = _nTextBoxComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_FIELD;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = propertyType;
	_nWidgets++;
	
	Label labelComponent;
	if (0 < label.length ()) {
	    labelComponent = new Label (label);
	} else {
	    labelComponent = new Label (name);
	}
	_panelPanelComponent.add (labelComponent);

	TextBox textBoxComponent 
	    = new TextBox ("", MAX_FIELD_CHARS, 1, 
			   TextBox.DISPLAY_ANY);
	_panelPanelComponent.add (textBoxComponent);
	_textBoxComponents[_nTextBoxComponents] = textBoxComponent;
	_nTextBoxComponents++;

	_panelHTMLLayout.br ();
    }

    private TextBox componentOfFieldWidget (int widgetId) {
	if (_widgetTypes[widgetId] != WIDGET_FIELD) {
	    error ("type must be field");
	}
	int componentIndex = _widgetComponentIndexs[widgetId];
	return _textBoxComponents[componentIndex];
    }



    /* choice widget */

    private void makeChoiceWidget (String name, String label) {
	int componentIndex = _nChoiceComponents;

	_widgetNames[_nWidgets] = name;
	_widgetTypes[_nWidgets] = WIDGET_CHOICE;
	_widgetComponentIndexs[_nWidgets] = componentIndex;
	_widgetPropertyTypes[_nWidgets] = PROPERTY_INTEGER;
	_nWidgets++;
	
	Label labelComponent;
	if (0 < label.length ()) {
	    labelComponent = new Label (label);
	} else {
	    labelComponent = new Label (name);
	}
	_panelPanelComponent.add (labelComponent);

	ListBox choiceComponent = new ListBox (ListBox.CHOICE);
	_panelPanelComponent.add (choiceComponent);
	_choiceComponents[_nChoiceComponents] 
	    = choiceComponent;
	_nChoiceComponents++;

	_panelHTMLLayout.br ();
    }

    private void addItemIntoChoiceWidget (String itemLabel) {
	ListBox currentChoiceComponent 
	    = _choiceComponents[_nChoiceComponents - 1];
	currentChoiceComponent.append (itemLabel);
    }

    private ListBox componentOfChoiceWidget (int widgetId) {
	if (_widgetTypes[widgetId] != WIDGET_CHOICE) {
	    error ("type must be choice");
	}
	int componentIndex = _widgetComponentIndexs[widgetId];
	return _choiceComponents[componentIndex];
    }



    /* application model I/O */

    private MediaImage _viewMediaImage;
    private Image _viewImage;

    // private String VIEW_IMAGE = "autogl_view_image.gif";
    private String VIEW_IMAGE = "autogl_view_image.jpeg";

    private String INFO_COMMAND_OUTPUT = "autogl_info_command_output_file";
    private String VALUE_COMMAND_OUTPUT = "autogl_value_command_output_file";
    private String POST_COMMAND_CGI = "autogl_post_command.cgi";
    
    private static final int VIEW_IMAGE_WIDTH = 240;
    private static final int VIEW_IMAGE_HEIGHT = 240;

    private void readImage () {
	if (_viewImage != null) {
	    _viewMediaImage.unuse ();
	    _viewImage.dispose ();
	    _viewMediaImage.dispose ();
	    _viewImage = null;
	    _viewMediaImage = null;
	}
	_viewMediaImage = MediaManager.getImage 
	    (getSourceURL () + VIEW_IMAGE);
	try {
	    _viewMediaImage.use ();
	    _viewImage = _viewMediaImage.getImage ();
	} catch (Exception ex) {
	    error (ex.toString ());
	    _viewImage = null;
	    _viewMediaImage = null;
	}
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
	    HttpConnection conn = (HttpConnection)Connector.open 
		(getSourceURL () + INFO_COMMAND_OUTPUT, 
		 Connector.READ);
	    conn.setRequestMethod (HttpConnection.GET);
	    conn.connect ();
	    initializeReader 
		(new InputStreamReader (conn.openInputStream ()));
	    parseInfoCommandOutputFile ();
	    terminateReader ();
	    conn.close ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
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
		ListBox component = componentOfToggleWidget (widgetId);
		if (value.equals ("true")) {
		    component.select (0);
		} else if (value.equals ("false")) {
		    component.deselect (0);
		} else {
		    error ("value must be true or false");
		    return;
		}
		
	    } else if (keyword.equals ("string")) {
		int widgetId = widgetIdOfName (name);
		TextBox component = componentOfFieldWidget (widgetId);
		component.setText (value);
		
	    } else if (keyword.equals ("integer")) {
		int widgetId = widgetIdOfName (name);
		if (_widgetTypes[widgetId] == WIDGET_FIELD) {
		    TextBox component = componentOfFieldWidget (widgetId);
		    component.setText (value);

		} else if (_widgetTypes[widgetId] == WIDGET_CHOICE) {
		    int integerValue = Integer.parseInt (value);
		    String itemLabel = readWord ();
		    ListBox component = componentOfChoiceWidget (widgetId);
		    component.select (integerValue);

		} else {
		    error ("integer must be field or choice");
		    return;
		}
		
	    } else if (keyword.equals ("real")) {
		int widgetId = widgetIdOfName (name);
		TextBox component = componentOfFieldWidget (widgetId);
		component.setText (value);
		
	    } else {
		error ("unsupported widget type");
		return;
	    }
	}
    }

    private void readModelValue () {
	try {
	    HttpConnection conn = (HttpConnection)Connector.open 
		(getSourceURL () + VALUE_COMMAND_OUTPUT, 
		 Connector.READ);
	    conn.setRequestMethod (HttpConnection.GET);
	    conn.connect ();
	    initializeReader 
		(new InputStreamReader (conn.openInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	    conn.close ();
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
		    ListBox component = 
			componentOfToggleWidget (iWidget);
		    if (component.isIndexSelected (0)) {
			writeSetCommand (writer, name, "true");
		    } else {
			writeSetCommand (writer, name, "false");
		    }
		}
		break;
	    case WIDGET_FIELD:
		{
		    TextBox component = 
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
		    ListBox component = 
			componentOfChoiceWidget (iWidget);
		    String value = "(NONE)";
		    for (int iItem = 0; 
			 iItem < component.getItemCount ();
			 iItem++) {
			if (component.isIndexSelected (iItem)) {
			    value = component.getItem (iItem);
			    break;
			}
		    }
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

    private void writeValueCommand (OutputStreamWriter writer) {
	try {
	    writer.write ("value");
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

    private void postCallCommand (String callbackName) {
	try {
	    HttpConnection conn = (HttpConnection)Connector.open 
		(getSourceURL () + POST_COMMAND_CGI, 
		 Connector.READ_WRITE);
	    conn.setRequestMethod (HttpConnection.POST);
	    conn.setRequestProperty ("Content-Type", "text/plain");
	    OutputStreamWriter writer = 
		new OutputStreamWriter (conn.openOutputStream ());
	    writeResizeCommand (writer);
	    writeSetCommandForEachWidget (writer);
	    writeCallCommand (writer, callbackName);
	    writeValueCommand (writer);
	    writer.close ();
	    conn.connect ();

	    initializeReader 
		(new InputStreamReader (conn.openInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	    conn.close ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void postDrawCommand () {
	try {
	    HttpConnection conn = (HttpConnection)Connector.open 
		(getSourceURL () + POST_COMMAND_CGI, 
		 Connector.READ_WRITE);
	    conn.setRequestMethod (HttpConnection.POST);
	    conn.setRequestProperty ("Content-Type", "text/plain");
	    OutputStreamWriter writer = 
		new OutputStreamWriter (conn.openOutputStream ());
	    writeResizeCommand (writer);
	    writeSetCommandForEachWidget (writer);
	    writeDrawCommand (writer);
	    writeValueCommand (writer);
	    writer.close ();
	    conn.connect ();

	    initializeReader 
		(new InputStreamReader (conn.openInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	    conn.close ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }

    private void postDragCommand 
	(int startX, int startY, int endX, int endY) {
	try {
	    HttpConnection conn = (HttpConnection)Connector.open 
		(getSourceURL () + POST_COMMAND_CGI, 
		 Connector.READ_WRITE);
	    conn.setRequestMethod (HttpConnection.POST);
	    conn.setRequestProperty ("Content-Type", "text/plain");
	    OutputStreamWriter writer = 
		new OutputStreamWriter (conn.openOutputStream ());
	    writeResizeCommand (writer);
	    writeSetCommandForEachWidget (writer);
	    writeDragCommand (writer, startX, startY, endX, endY);
	    writeDrawCommand (writer);
	    writeValueCommand (writer);
	    writer.close ();
	    conn.connect ();

	    initializeReader 
		(new InputStreamReader (conn.openInputStream ()));
	    parseValueCommandOutputFile ();
	    terminateReader ();
	    conn.close ();
	} catch (IOException ex) {
	    error (ex.toString ());
	}
    }



    /* public methods */

    public Image viewImage () { 
	return _viewImage;
    }
 
    public void switchToPanel () {
	Display.setCurrent (_panelPanelComponent);
    }

    public void drawView () {
	postDrawCommand ();
	readImage ();
    }

    public void dragView 
	(int startX, int startY, int endX, int endY) {
	postDragCommand (startX, startY, endX, endY);
	readImage ();
    }

    /* overwride */
    public void componentAction (Component source, int type, int param) {
	for (int iComponent = 0; 
	     iComponent < _nButtonComponents; 
	     iComponent++) {
	    if (source == _buttonComponents[iComponent]
		&& type == BUTTON_PRESSED) {
		int widgetId = _buttonWidgetIds[iComponent];
		String name = _widgetNames[widgetId];
		
		postCallCommand (name);
		readImage ();
	    }
	}
    }

    /* overwride */
    public void softKeyPressed (int key) {
    }

    /* overwride */
    public void softKeyReleased (int key) {
	if (key == Frame.SOFT_KEY_2) {
	    switchToView ();
	}
    }

    /* overwride */
    public void start () {
	clearWidget ();
	makePanel ();
	makeView ();
	switchToPanel ();
    }
}

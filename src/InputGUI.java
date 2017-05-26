import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.TreeMap;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class InputGUI extends JFrame implements ActionListener
{
	final int[] dimension = {1000, 800}; // fixed GUI dimension
	TreeMap<String, String> parameterValues;
	
	JPanel[] inputPanels = new JPanel[9];
	JPanel[] commentPanels = new JPanel[inputPanels.length];
	
	// file chooser button
	JFileChooser fc = new JFileChooser();
	JButton candidateButton;
	JButton mzXMLButton;
	JButton dbButton;
	JButton dbReverseButton;
	JButton outputNameButton;
	JButton paramButton;
	
	// text input
	JTextField candidatePath = new JTextField(50);
	JTextField mzXMLPath = new JTextField(50);
	JTextField dbPath = new JTextField(20);
	JTextField dbReversePath = new JTextField(20);
	JTextField outputName = new JTextField(50);
	JTextField paramPath = new JTextField(50);
	
	JTextField halfRTWidth = new JTextField(10);
	JTextField minGMscore = new JTextField(10);
	
	JTextField minIonCov = new JTextField(10);
	JTextField minDelta = new JTextField(10);
	JTextField maxFDR = new JTextField(10);
	
	// combo box
	String[] numCPUchoices;
	JComboBox<String> numCPU;
	
	// check box
	JCheckBox targetList = new JCheckBox("Output target list only.");
	
	// control
	JButton submit = new JButton("Process");
	JButton clear = new JButton("Clear");
	
	public InputGUI()
	{
		super("ZXMiner 1.0.3 - April 2017");
		
		setPreferredSize(new Dimension(dimension[0], dimension[1]));
		setLayout(new FlowLayout());
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		setup();
		
		Container contentPane = getContentPane();
		
		for (int i = 0; i < inputPanels.length; i++)
		{
			contentPane.add(inputPanels[i]);
			contentPane.add(commentPanels[i]);
		}
		
		pack();
        setVisible(true);
	}
	
	public void setup()
	{
		for (int i = 0; i < inputPanels.length; i++)
		{
			inputPanels[i] = new JPanel(new FlowLayout(FlowLayout.LEFT));
			inputPanels[i].setPreferredSize(new Dimension(dimension[0] - 100, 35));
			
			commentPanels[i] = new JPanel(new FlowLayout(FlowLayout.LEFT));
			commentPanels[i].setPreferredSize(new Dimension(dimension[0] - 100, 35));
		}
		
		// candidate path
		inputPanels[0].add(new JLabel("> candidate precursor file: "));
		candidatePath.setEditable(false);
		inputPanels[0].add(candidatePath);
		candidateButton = new JButton("Browse...");
		candidateButton.addActionListener(this);
		inputPanels[0].add(candidateButton);
		commentPanels[0].add(new JLabel("   'Select path to precursor file. Default location is \".\\input\\\".'"));
		
		// mzXML path
		inputPanels[1].add(new JLabel("> mzXML folder: "));
		mzXMLPath.setEditable(false);
		inputPanels[1].add(mzXMLPath);
		mzXMLButton = new JButton("Browse...");
		mzXMLButton.addActionListener(this);
		inputPanels[1].add(mzXMLButton);
		commentPanels[1].add(new JLabel("   'Select directory path to mzXML files.'"));
		
		// fasta path
		inputPanels[2].add(new JLabel("> forward fasta file: "));
		dbPath.setEditable(false);
		inputPanels[2].add(dbPath);
		dbButton = new JButton("Browse...");
		dbButton.addActionListener(this);
		inputPanels[2].add(dbButton);
		inputPanels[2].add(new JLabel("> reverse fasta file: "));
		inputPanels[2].add(dbReversePath);
		dbReverseButton = new JButton("Browse...");
		dbReverseButton.addActionListener(this);
		inputPanels[2].add(dbReverseButton);
		commentPanels[2].add(new JLabel("   'Select fasta file(s). Default location is \".\\sequence\\\". Use reversed DB only for targeted stage.'"));
		
		// output file name
		inputPanels[3].add(new JLabel("> output file name: "));
		inputPanels[3].add(outputName);
		outputNameButton = new JButton("Browse...");
		outputNameButton.addActionListener(this);
		inputPanels[3].add(outputNameButton);
		commentPanels[3].add(new JLabel("   'Manually enter output file name or select an existing \"_allScores.out\" file in \".\\output\\\" to re-process the result without performing full re-analysis.'"));
		
		// candidate path
		inputPanels[4].add(new JLabel("> configuration file: "));
		paramPath.setEditable(false);
		inputPanels[4].add(paramPath);
		paramButton = new JButton("Browse...");
		paramButton.addActionListener(this);
		inputPanels[4].add(paramButton);
		commentPanels[4].add(new JLabel("   'Select configuration file. Default location is \".\\config\\\".'"));
		
		// target list
		inputPanels[5].add(new JLabel("> target list parameters: "));
		inputPanels[5].add(targetList);
		inputPanels[5].add(new JLabel("> half retention time width (min): "));
		inputPanels[5].add(halfRTWidth);
		inputPanels[5].add(new JLabel("> minimum GM score: "));
		minGMscore.setText("0.5");
		inputPanels[5].add(minGMscore);
		commentPanels[5].add(new JLabel("   'Specify parameters for outputing target list (for high-res MS/MS). Using 4 min for 4hr gradient or 2 min for 85m gradient is recommended'"));
		
		// target list
		inputPanels[6].add(new JLabel("> minimum ion coverage score: "));
		minIonCov.setText("0.2");
		inputPanels[6].add(minIonCov);
		inputPanels[6].add(new JLabel("> minimum delta score: "));
		minDelta.setText("0.1");
		inputPanels[6].add(minDelta);
		inputPanels[6].add(new JLabel("> maximum FDR (decimal): "));
		maxFDR.setText("0.01");
		inputPanels[6].add(maxFDR);
		commentPanels[6].add(new JLabel("   'Parameters for evaluating high-resolution MS/MS data. Not used when only creating target list.'"));
		
		// num CPU
		numCPUchoices = new String[Runtime.getRuntime().availableProcessors()];
		
		for (int i = 0; i < numCPUchoices.length; i++)
			numCPUchoices[i] = (i + 1) + "";
		
		numCPU = new JComboBox<String>(numCPUchoices);
		
		inputPanels[7].add(new JLabel("> number of CPU cores: "));
		inputPanels[7].add(numCPU);
		commentPanels[7].add(new JLabel("   'Select number of CPU cores to allocate to this project.'"));
		
		// control
		inputPanels[8].add(submit);
		inputPanels[8].add(clear);
		submit.addActionListener(this);
		clear.addActionListener(this);
	}

	// actions
	public void actionPerformed(ActionEvent a)
	{
		if (a.getSource().equals(candidateButton))
		{
			fc.setCurrentDirectory(new File(".\\input\\"));
			fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			int returnVal = fc.showOpenDialog(this);

        	if (returnVal == JFileChooser.APPROVE_OPTION)
        	{
        		File file = fc.getSelectedFile();
        		candidatePath.setText(".\\input\\" + file.getName());
        	}

			return;
		}
		
		if (a.getSource().equals(mzXMLButton))
		{
			fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int returnVal = fc.showOpenDialog(this);

        	if (returnVal == JFileChooser.APPROVE_OPTION)
        	{
        		File file = fc.getSelectedFile();
        		mzXMLPath.setText(file.getPath());
        	}

			return;
		}
		
		if (a.getSource().equals(dbButton))
		{
			fc.setCurrentDirectory(new File(".\\sequence\\"));
			fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			int returnVal = fc.showOpenDialog(this);

        	if (returnVal == JFileChooser.APPROVE_OPTION)
        	{
        		File file = fc.getSelectedFile();
        		dbPath.setText(file.getName());
        	}

			return;
		}
		
		if (a.getSource().equals(dbReverseButton))
		{
			fc.setCurrentDirectory(new File(".\\sequence\\"));
			fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			int returnVal = fc.showOpenDialog(this);

        	if (returnVal == JFileChooser.APPROVE_OPTION)
        	{
        		File file = fc.getSelectedFile();
        		dbReversePath.setText(file.getName());
        	}

			return;
		}
		
		if (a.getSource().equals(outputNameButton))
		{
			fc.setCurrentDirectory(new File(".\\output\\"));
			fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			int returnVal = fc.showOpenDialog(this);

        	if (returnVal == JFileChooser.APPROVE_OPTION)
        	{
        		File file = fc.getSelectedFile();
        		outputName.setText(file.getPath());
        	}

			return;
		}
		
		if (a.getSource().equals(paramButton))
		{
			fc.setCurrentDirectory(new File(".\\config\\"));
			fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			int returnVal = fc.showOpenDialog(this);

        	if (returnVal == JFileChooser.APPROVE_OPTION)
        	{
        		File file = fc.getSelectedFile();
        		paramPath.setText(file.getPath());
        	}

			return;
		}
		
		if (a.getSource().equals(clear))
		{
			candidatePath.setText("");
			mzXMLPath.setText("");
			dbPath.setText("");
			dbReversePath.setText("");
			outputName.setText("");
			paramPath.setText("");

			return;
		}
		
		if (a.getSource().equals(submit))
		{
			parameterValues = new TreeMap<String, String>();
			parameterValues.put("candidateFilePath", candidatePath.getText());
			parameterValues.put("mzXMLFilePath", mzXMLPath.getText());
			parameterValues.put("forwardFastaFileName", dbPath.getText());
			parameterValues.put("decoyFastaFileName", dbReversePath.getText());
			parameterValues.put("outputFileName", outputName.getText());
			parameterValues.put("paramFilePath", paramPath.getText());
			
			if (targetList.isSelected())
				parameterValues.put("outputTargetList", "true");
			else
				parameterValues.put("outputTargetList", "false");
			
			parameterValues.put("targetListRTHalfWindow", halfRTWidth.getText());
			parameterValues.put("minScoreTargetList", minGMscore.getText());
			
			parameterValues.put("minIonCov", minIonCov.getText());
			parameterValues.put("minDeltaScore", minDelta.getText());
			parameterValues.put("maximumFDR", maxFDR.getText());
			
			parameterValues.put("numCPU", (String) numCPU.getSelectedItem());
			
			new XlinkMiner(new ParamStruct(parameterValues));
			
			return;
		}
	}
	
	public static void main(String[] args)
	{
		new InputGUI();
		
		// ParamStruct param = new ParamStruct();
		// CandidatePicker.pickPeptide(param);
	}
}

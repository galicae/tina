package util;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ContainerEvent;
import java.awt.event.ContainerListener;
import java.awt.event.HierarchyBoundsListener;
import java.awt.event.HierarchyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.HashMap;


import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolSimpleViewer;
import org.openscience.jmol.app.Jmol;
import org.openscience.jmol.app.jmolpanel.JmolPanel;

import bioinfo.proteins.PDBEntry;

public class JMolView {

	private String title = "Java induced JMol-Frame";
	private JmolPanel jmolPanel;
	private JFrame frame;
	private JmolSimpleViewer viewer;

	

	public JMolView() {
		frame = new JFrame();
		frame.setDefaultCloseOperation(frame.EXIT_ON_CLOSE);
		
		Container contentPane = frame.getContentPane();
		jmolPanel = new JmolPanel();

		jmolPanel.setPreferredSize(new Dimension(400, 400));
		frame.setResizable(false);
		contentPane.add(jmolPanel);

		frame.pack();
		frame.setVisible(true);
	}
	
	public void executeCmd(String cmd){
		jmolPanel.executeCmd(cmd);
	}
	
	public void addPDBEntry(PDBEntry entry, String color){
		jmolPanel.addPDBEntry(entry, color);
	}

	static class JmolPanel extends JPanel {
		
		private static String content= "";
		private static int nextContent = 1;
		private static HashMap<Integer,String> color = new HashMap<Integer,String>();

		private static final long serialVersionUID = 1L;
		private static JmolSimpleViewer viewer;
		private static JmolAdapter adapter;

		public JmolPanel() {
			adapter = new SmarterJmolAdapter();
			viewer = JmolSimpleViewer.allocateSimpleViewer(this, adapter);
		}

		public JmolSimpleViewer getViewer() {
			return viewer;
		}

		public void executeCmd(String rasmolScript) {
			viewer.evalString(rasmolScript);
		}

		final Dimension currentSize = new Dimension();
		final Rectangle rectClip = new Rectangle();

		public void paint(Graphics g) {
			getSize(currentSize);
			//g.getClipBounds(rectClip);
			viewer.renderScreenImage(g, (int)currentSize.getHeight(), (int)currentSize.getWidth());
		}
		
		public void addPDBEntry(PDBEntry entry, String jmolColor){
			color.put(nextContent, jmolColor);
			content += "MODEL     "+String.format("%4d",nextContent++)+"\n";
			content += entry.getAtomSectionAsString();
			content += "ENDMDL\n";
			//System.out.println(content+"\n\n");
			
			viewer.openStringInline(content);
			executeCmd("cartoon ONLY");
			for(int i = 1; i < nextContent; i++){
				//System.out.println(i+" "+color.get(i));
				executeCmd("model "+i);
				executeCmd("select visible");
				executeCmd("color "+color.get(i));
			}
			executeCmd("model 0");
		}
	}

}

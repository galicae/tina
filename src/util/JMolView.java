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


import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolSimpleViewer;
import org.openscience.jmol.app.Jmol;
import org.openscience.jmol.app.jmolpanel.JmolPanel;

public class JMolView {

	private String title = "Java induced JMol-Frame";
	private JmolPanel jmolPanel;
	private JFrame frame;
	private JmolSimpleViewer viewer;
	

	public JMolView(String inlinePdb) {
		frame = new JFrame();
		frame.setDefaultCloseOperation(frame.EXIT_ON_CLOSE);
		
		Container contentPane = frame.getContentPane();
		jmolPanel = new JmolPanel(inlinePdb);

		jmolPanel.setPreferredSize(new Dimension(400, 400));
		frame.setResizable(false);
		contentPane.add(jmolPanel);

		frame.pack();
		frame.setVisible(true);
	}

	static class JmolPanel extends JPanel {

		private static final long serialVersionUID = 1L;
		private static JmolSimpleViewer viewer;
		private static JmolAdapter adapter;

		public JmolPanel(String inlinePdb) {
			adapter = new SmarterJmolAdapter();
			viewer = JmolSimpleViewer.allocateSimpleViewer(this, adapter);
			viewer.openStringInline(inlinePdb);
			executeCmd("cartoon ONLY");
			executeCmd("model 1");
			executeCmd("select visible");
			executeCmd("color red");
			executeCmd("model 2");
			executeCmd("select visible");
			executeCmd("color green");
			executeCmd("model 0");
			
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
	}

	public static void main(String[] args) {
		String inlinePdb = "ATOM     31  CB  VAL A   4      14.588  -8.758  20.143";
			JMolView jmol = new JMolView(inlinePdb);			
			
	}

}

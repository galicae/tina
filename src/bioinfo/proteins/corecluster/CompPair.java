package bioinfo.proteins.corecluster;

import java.io.Serializable;

public class CompPair<T extends Comparable<T>> implements Serializable, Comparable<CompPair<T>> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private T x, y;

	public CompPair(T x, T y) {
		this.x = x;
		this.y = y;
	}

	public T getX() {
		return x;
	}

	public void setX(T x) {
		this.x = x;
	}

	public T getY() {
		return y;
	}

	public void setY(T y) {
		this.y = y;
	}

	@Override
	public String toString() {
		return "(" + x + " | " + y + ")";
	}

	@Override
	public int compareTo(CompPair<T> o) {
		if (this.x.compareTo(o.getY()) > 0) {
			return 1;
		} else if (this.y.compareTo(o.getX()) < 0) {
			return -1;
		} else {
			return 0;
		}
	}

}

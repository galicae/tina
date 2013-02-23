package bioinfo.proteins.corecluster;

import java.io.Serializable;

public class Pair<T> implements Serializable {

	private static final long serialVersionUID = 3964318084769219061L;
	private T x, y;

	public Pair() {
	}

	public Pair(T x, T y) {
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

}

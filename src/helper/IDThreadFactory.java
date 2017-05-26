import java.util.concurrent.ThreadFactory;

// define thread factory which provide predictable 'prefix-n' as thread name
// n = {0, 1, ..., numCPU - 1}
public class IDThreadFactory implements ThreadFactory
{
	private int counter = 0;
	private String prefix = "";

 	public IDThreadFactory(String prefix)
 	{ this.prefix = prefix; }

 	public Thread newThread(Runnable r)
 	{ return new Thread(r, prefix + "-" + counter++); }
}
/*
 * Simplest possible test - just execute a couple of commands and 
 * join them together in a pipeline
 */
hello = {
	exec "echo hello > test.txt"
  exec "sleep 6"
}

world = {
	exec "echo world > test.world.txt"
  exec "sleep 6"
}

Bpipe.run {
	hello + world 
}

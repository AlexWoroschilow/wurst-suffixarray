sub move($ $ $ $){
	my($from, $to, $via, $count) = @_;
	if($count > 1){
		move($from, $via, $to, $count-1);
		move($from, $to, $via, 1);
		move($via, $to, $from, $count-1);
	}
	else{
		print "move disc from $from to $to\n";
	}
}
move($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

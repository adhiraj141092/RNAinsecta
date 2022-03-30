#!/usr/bin/perl -w
use warnings;
use strict;
my $inFile="&STDIN";
my $outFile="&STDOUT";

open (inFile, "<$inFile") or die( "Cannot open input file $inFile: $!" );
open (outFile, ">$outFile") or die ("Cannot open output file $outFile: $!");

my @weight = ();
get_weight_table(\@weight);

$method = 1;
$reverse_flag = 0;
$num_main=0;
$num_total=0;

print outFile "A... A..( A.(. A.(( A(.. A(.( A((. A((( G... G..( G.(. G.(( G(.. G(.( G((. G((( C... C..( C.(. C.(( C(.. C(.( C((. C((( U... U..( U.(. U.(( U(.. U(.( U((. U(((\n";
while(<inFile>){
	if($_ =~ />(.+)/){	
		$num_total++;
		$iden = $1;
		next;
	}
	if(uc($_) =~ /([ACUG]+)/){	
		$seq_letter = $1;
		print "$seq_letter\n\n";
		next;
	}
	if($_ =~ /(.*?)\s+/){
		$symble = $1;	
		@struct_symble = divide_structure($symble);
		@main_struct = get_main_struct(\@struct_symble);
		# print "$iden\n\n";
		# print "$symble\n\n";
		# print "@struct_symble\n\n";
		# print "@main_struct\n\n";
		if(scalar(@main_struct) > 1){
			$num_main++;
			$left_arm = get_match_all_body($main_struct[0]);
			$right_arm = get_match_all_body($main_struct[1]);
			$left_posi = get_position_in_symble($left_arm, $symble);
			$right_posi = get_position_in_symble($right_arm, $symble);
			$left_raw_letter = get_seq_in_raw($seq_letter, $left_posi,length($left_arm));
			$right_raw_letter = get_seq_in_raw($seq_letter, $right_posi, length($right_arm));

			if($reverse_flag == 0){		
				@coding_table = translate_to_coding($left_raw_letter, $left_arm, $right_raw_letter, $right_arm);
			}
			else{
				$reverse_right_arm = reverse($right_arm);
				$reverse_right_raw_letter = reverse($right_raw_letter);
				@coding_table = translate_to_coding($left_raw_letter, $left_arm, $reverse_right_raw_letter, $reverse_right_arm);
			}
			
			$len_coding_table = @coding_table;

			if($method == 1){
				unite_probability_(\@coding_table);					
			}
			elsif($method == 2){
				bool_(\@coding_table);				
			}
			elsif($method == 3){
				product_weight_(\@coding_table);
			}
			else{
				print "ERROR: not choose method\n";
			}			

			unite_probability_(\@coding_table);
			for($i=0; $i<$len_coding_table; $i++){
				print outFile "$coding_table[$i] ";	
			}
			print outFile "\n";
		}		
	}

}

print "$num_total queries, triplet coding $num_main queries!\n";
#######################################
# sub fuctions
sub get_weight_table{
	local($table) = @_;
	for(my $i=0; $i<32; $i+=8){
		$$table[$i+0] = 0.15 / 8;
		$$table[$i+1] = 0.05 / 8;
		$$table[$i+2] = 0.25 / 8;
		$$table[$i+3] = 0.15 / 8;
		$$table[$i+4] = 0.05 / 8;
		$$table[$i+5] = 0.1 / 8;
		$$table[$i+6] = 0.15 / 8;
		$$table[$i+7] = 0.1 / 8;		
		
	}
	
	return @$table;
}

sub product_weight_{
	local($table) = @_;
	for my $i (0..@$table-1) {
		$$table[$i] *= $weight[$i];			
	}
}

sub bool_{
	local($table) = @_;
	for my $i (0..@$table-1) {
		$$table[$i] = 1 if($$table[$i] != 0);	
	}
}

sub unite_probability_{
	local($table) = @_;
	
	my $sum = 0;
	for my $i (0..@$table-1) {
		$sum += $$table[$i];			
	}

	for my $i (0..@$table-1) {
		$$table[$i] /= $sum;	
	}
}

sub get_char_value{
	local($char) = @_;

	my $l_char = lc($char);
	if($l_char eq 'a'){
		return 0;	
	}
	elsif($l_char eq 'g'){
		return 1;	
	}
	elsif($l_char eq 'c'){
		return 2;	
	}
	elsif($l_char eq 'u' || $l_char eq 't'){
		return 3;	
	}
	elsif($l_char eq 's'){
		return 4;	
	}
	else{
		print "ERROR: has not A G C U or S(deleted stem-loop)\n";
		return -1;
	}
}

sub get_near_value{
	local($line) = @_;
	
	my $len_line = length $line;
	my $value = 0;
	for(my $i=$len_line-1; $i>-1; $i--){
		$value += 2**($len_line - 1 - $i) * get_dot_brackle_value(substr($line, $i, 1));	
	}
	
	return $value;
}

sub get_dot_brackle_value{
	local($char) = @_;
	if($char eq '.'){
		return 0;	
	}
	elsif($char =~ m/[\(\)]/) {
		return 1;	
	}
	else{
		print "ERROR: has neither . or ( or )\n";
		return -1;
	}
}

sub translate_to_coding{
	local($left_seq, $left_arm, $right_seq, $right_arm) = @_;
	local($len_left_arm, $len_right_arm, $i, $char, $value_char, $value_near);
	local(@table);
	
	foreach(0..31){
		push(@table, 0);	
	}
	
	$len_left_arm = length $left_arm;
	for($i=0; $i<$len_left_arm; $i++){
		$value_char = get_char_value(substr($left_seq, $i, 1));
		return @table if($value_char == -1);	
		
		
		if($i == 0){
			$value_near = get_near_value(".".substr($left_arm, $i, 2));
		}	
		elsif($i == ($len_left_arm - 1)){
			$value_near = get_near_value(substr($left_arm, $i-1, 2).".");
		}
		else{
			$value_near = get_near_value(substr($left_arm, $i-1, 3));
		}
		
		$table[$value_char * 8 + $value_near]++;
	}

	$len_right_arm = length $right_arm;	
	for($i=0; $i<$len_right_arm; $i++){
		$value_char = get_char_value(substr($right_seq, $i, 1));
	#	print outFile "$i	char is $char, value is $value_char	";
	#	if($value_char == -1){
	#		return @table;	
	#	}
		
		
		if($i == 0){
			$value_near = get_near_value(".".substr($right_arm, $i, 2));
		}	
		elsif($i == ($len_right_arm - 1)){
			$value_near = get_near_value(substr($right_arm, $i-1, 2).".");
		}
		else{
			$value_near = get_near_value(substr($right_arm, $i-1, 3));
		}
		$table[$value_char * 8 + $value_near]++;
	}
	
	return @table;
}

sub get_seq_in_raw{
	local($seq, $posi, $len_arm) = @_;
	return substr($seq, $posi, $len_arm);
}

sub get_position_in_symble{
	local($str, $seq) = @_;	
	return (index($seq, $str));
}

sub get_match_all_body{
	local($seq) = @_;
	local($len_seq, $i, $char, $ext_str, $len_ext_str, $body);
	
	$len_seq = length $seq;
	for($i=0; $i<$len_seq; $i++){
		if(substr($seq, $i, 1) ne "\."){
			$ext_str = substr($seq, $i);
			last;	
		}					
	}
	$len_ext_str = length $ext_str;
	$body = $ext_str;
	for($i=($len_ext_str-1); $i>-1; $i--){
		(substr($ext_str, $i, 1) eq "\.") ? chop($body) : last;
	}
	
#	print outFile "$body\n";
	
	return $body;			
}

sub get_left_bracket_number{
	local($seq) = @_;
	return ($seq =~ tr/'('//);
}

sub get_right_bracket_number{
	local($seq) = @_;
	return ($seq =~ tr/')'//);
}

sub get_main_struct{
	local($struct) = @_;
	local($len_struct, $i, $len_left_bracket); 
	local(@left_bracket, @right_bracket);
	
	$len_struct = @$struct;
	for($i=0; $i<$len_struct; $i+=2){
		push(@left_bracket, get_left_bracket_number($$struct[$i]));
		push(@right_bracket, get_right_bracket_number($$struct[$i+1]));
	#	print "left is $left_bracket_number, right is $right_bracket_number\n";
	}
	
	$same=0;
	$len_left_bracket = @left_bracket;
	for($i=0; $i<$len_left_bracket; $i++){
		$same++ if($left_bracket[$i] == $right_bracket[$i]);
	}
	
	return ($same >= 1) ? @$struct : ();
}

sub divide_structure{
	local($symble) = @_;
	local($bracket_flag, $len_symble, $i, $char, @posi, $len_posi);
	local(@struct);
	
	push(@posi, 0);
	$bracket_flag = '(';
	$len_symble = length $symble;
	for($i=0; $i<$len_symble; $i++){
		$char = substr($symble, $i, 1);
		if($char ne '.' && $char ne $bracket_flag){
			$bracket_flag = $char;
			push(@posi, $i);
		}
	}
	
	$len_posi = @posi;
	for($i=0; $i<($len_posi-1); $i++){
		push(@struct, substr($symble, $posi[$i], $posi[$i+1] - $posi[$i]));
	}

	push(@struct, substr($symble, $posi[$len_posi-1]));
	return @struct;
}

###############################################
close inFile or die "can't close the input file : $!";
close outFile or die "can't close the output file : $!"; 

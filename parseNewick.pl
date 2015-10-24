#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $leaf = qr/([A-Z]+)/;
my $node = qr/(\(.+\))/;
my $leaf_leaf = qr/^\($leaf,$leaf\)$/;
my $node_leaf = qr/^\($node,$leaf\)$/;
my $leaf_node = qr/^\($leaf,$node\)$/;
my $node_node = qr/^\($node,$node\)$/;

my $nt = "((A,(B,((C,(D,E)),(F,G)))),H)";
my $tree = newick ($nt,undef);
print Dumper $tree;

sub newick {
  my ($t,$tree) = @_;
  if ($t =~ /$leaf_node/){
    my ($leaf,$node) = ($1,$2);
    ins_Leaf_Node ($tree,$leaf, newick ($node,$tree));
  } elsif ($t =~ /$node_leaf/){
    my ($node,$leaf) = ($1,$2);
    ins_Node_Leaf ($tree,newick ($node,$tree),$leaf);
  } elsif ($t =~ /$node_node/){
    my ($node1,$node2) = ($1,$2);
    ins_Node_Node ($tree,newick ($node1,$tree), newick ($node2,$tree));
  } elsif ($t =~ /$leaf_leaf/){
    my ($leaf1,$leaf2) = ($1,$2);
    ins_Leaf_Leaf ($tree,$leaf1,$leaf2);
  } else {
    die "Unrecognized tree branch $t\n";
  }
}

sub create_leaf {
  my ($value) = @_;
  return {
	  'VALUE' => $value,
	 };
}

sub ins_Node_Node {
  my ($tree,$left,$right) = @_;
  $tree->{LEFT} = $left;
  $tree->{RIGHT} = $right;
  return $tree;
}

sub ins_Node_Leaf {
  my ($tree,$left,$right) = @_;
  $tree->{LEFT} = $left;
  $tree->{RIGHT} = create_leaf ($right);
  return $tree;
}

sub ins_Leaf_Node {
  my ($tree,$left,$right) = @_;
  $tree->{LEFT} = create_leaf ($left);
  $tree->{RIGHT} = $right;
  return $tree;
}

sub ins_Leaf_Leaf {
  my ($tree,$left,$right) = @_;
  my $nodeR = create_leaf ($right);
  my $nodeL = create_leaf ($left);
  $tree->{LEFT}  = $nodeL;
  $tree->{RIGHT} = $nodeR;
  return $tree;
}



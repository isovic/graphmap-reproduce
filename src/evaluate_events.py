#! /usr/bin/python

import os;
import sys;
import math;
import random;

class StructEvent:
	def __init__(self, size=0, type='', position=0):
		self.size = size;
		self.type = type;
		self.position = position;
		self.position_final = position;
		self.open_cov = 0;
		self.close_cov = 0;
		self.max_cov = 0;
		self.max_event_length = 0;
		self.comment = '';
		self.evaluation = '(not evaluated)';

	def verbose(self):
		sys.stderr.write(self.verbose_to_string());

	def verbose_to_string(self):
		event_type = 'insertion_in_read' if (self.type == 'D' or self.type == 'insertion_in_read') else 'deletion_in_read';
		return 'Event type: %s, position: %7d, end: %7d, size: %4d, max_event_length = %4d, result: %s, comment: %s' % (event_type, self.position, (self.position + self.size), self.size, self.max_event_length, self.evaluation, self.comment);

	def csv_line(self):
		# ret = 'type: %s\tstart: %d\tend: %d\topen_cov: %d\tclose_cov: %d\tmax_cov: %d' % (self.type, self.start, self.end, self.open_coverage, self.close_coverage, self.max_coverage);
		event_type = 'insertion_in_read' if (self.type == 'D' or self.type == 'insertion_in_read') else 'deletion_in_read';
		ret = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s' % (event_type, self.position, (self.position + self.size), self.size, self.open_cov, self.close_cov, self.max_cov, self.max_event_length, self.comment);
		return ret;

	def csv_header(self):
		ret = 'type\tstart\tend\tlength\topen_cov\tclose_cov\tmax_cov\tmax_event_length\tcomment';
		return ret;

def check_event_near_blacklist(event, blacklist_events, allowed_distance):
	# print event.verbose_to_string();
	for blacklist_event in blacklist_events:
		if (blacklist_event.comment != 'real'):
			# print blacklist_event.verbose_to_string();

			event_start = event.position;
			event_end = event.position + event.size;
			blacklist_event_start = blacklist_event.position;
			# blacklist_event_end = (blacklist_event.position + blacklist_event.size) if (blacklist_event.type == 'I') else (blacklist_event.position + 1);
			blacklist_event_end = (blacklist_event.position + blacklist_event.size);
			# print 'Blacklisted event: ', blacklist_event.verbose_to_string();
			# print blacklist_event_start;
			# print blacklist_event_end;

			if ((event_start < blacklist_event_start and event_end < blacklist_event_start) or
				(event_start > blacklist_event_end and event_end > blacklist_event_end)):

				### If the events do not overlap, check the distances.
				dists = [abs(event_start - blacklist_event_start),
						 abs(event_start - blacklist_event_end),
						 abs(event_end - blacklist_event_start),
						 abs(event_end - blacklist_event_end)];
				# print dists;
				for dist in dists:
					if (dist < allowed_distance):
						return True;
			else:
				### If they overlap, just return True.
				return True;

	# if (event.position == 1096291):
	# 	print 'Tu sam!';
	# 	print allowed_distance;
	# 	exit(1);

	return False;

def evaluate_event(event, truth_events, allowed_percent_variation):
	event_distances = [];

	allowed_distance = event.size * allowed_percent_variation;

	for truth_event in truth_events:
		event_start = event.position;
		event_end = event.position + event.size;
		truth_event_start = truth_event.position;
		truth_event_end = (truth_event.position + truth_event.size) if (truth_event.type == 'I') else (truth_event.position + 1);
#####		allowed_distance = truth_event.size * allowed_percent_variation;

		if (check_event_near_blacklist(event, truth_events, allowed_distance) == True):
			event.evaluation = TN;
			return TN;
		else:
			distance = abs(truth_event_start - event_start);
			event_distances.append([distance, truth_event]);

	sorted_event_distances = sorted(event_distances, key=lambda x: x[0]);

	### Take only the truth event with the smallest distance.
	### best_truth[0] is the distance from the tested event to the truth event.
	### best_truth[1] is a reference to the truth event.
	best_truth = sorted_event_distances[0];
	truth_size = best_truth[1].size;
	size_difference = abs(event.max_event_length - truth_size);
#####	allowed_distance = best_truth[1].size * allowed_percent_variation;
	# sys.stderr.write(best_truth[1].verbose_to_string() + '\n');

	### Check the starting position of the event, as well as its length.
	if (best_truth[1].type == event.type and best_truth[0] < allowed_distance and size_difference < math.ceil(allowed_percent_variation * truth_size)):
		event.evaluation = TP;
		best_truth_type = 'insertion_in_read' if (best_truth[1].type == 'D' or best_truth[1].type == 'insertion_in_read') else 'deletion_in_read';
		event.comment += ' nearest_truth=[dist:%d size:%d pos:%d type:%s]' % (best_truth[0], best_truth[1].size, best_truth[1].position, best_truth_type);
		return TP;

	# print '-> %s\t%s' % (best_truth[0], best_truth[1].verbose_to_string());

	event.evaluation = FP;

	# event.evaluation += ' nearest_truths=[%s]' % (';'.join([str(truth_event[1].position) for truth_event in sorted_event_distances]));
	# event.comment += ' nearest_truth=[dist:%d type:%s pos:%d size:%d]' % (best_truth[0], best_truth[1].type, best_truth[1].position, best_truth[1].size);
	best_truth_type = 'insertion_in_read' if (best_truth[1].type == 'D' or best_truth[1].type == 'insertion_in_read') else 'deletion_in_read';
	event.comment += ' nearest_truth=[dist:%d type:%s pos:%d size:%d]' % (best_truth[0], best_truth_type, best_truth[1].position, best_truth[1].size);

	if (best_truth[0] >= allowed_distance):
		event.comment += ' dist=%d/%d' % (best_truth[0], allowed_distance);
	if (size_difference >= allowed_percent_variation * truth_size):
		event.comment += ' size=%d/%d' % (size_difference, math.ceil(allowed_percent_variation * truth_size));
	if (best_truth[1].type != event.type):
		event.comment += ' wrong_type';

	return FP;

def evaluate(events, events_truth, fp_out=sys.stdout):
	sorted_events = sorted(events, key=lambda x: x.position);
	sorted_events_truth = sorted(events_truth, key=lambda x: x.position);

	i = 0;
	for event in sorted_events:
		i += 1;
		fp_out.write('[%d] %s\t%s\n' % (i, evaluate_event(event, events_truth, 0.25), event.verbose_to_string()));
		# print '\n';



def load_events(events_path):
	events = [];

	num_lines = 0;
	fp = open(events_path);
	for line in fp:
		num_lines += 1;
		line = line.strip();
		if (len(line) == 0):
			continue;

		if (num_lines == 1):
			labels = line;
			continue;

		# print line;

		values = line.split('\t');
		# event_type = 'I' if (values[0] == 'deletion_in_read') else 'D';
		event_type = 'insertion_in_read' if (values[0] == 'D' or values[0] == 'insertion_in_read') else 'deletion_in_read';
		start = int(values[1]);
		end = int(values[2]);
		length = int(values[3]);
		open_cov = int(values[4]);
		close_cov = int(values[5]);
		max_cov = int(values[6]);
		max_event_length = int(values[7]);
		comment = values[8] if (len(values) > 8) else '';

		event = StructEvent(length, event_type, start);
		event.open_cov = open_cov;
		event.close_cov = close_cov;
		event.max_cov = max_cov;
		event.max_event_length = max_event_length;
		event.comment = comment;

		events.append(event);

	fp.close();

	### Chain events that got fragmented
	sorted_events = sorted(events, key=lambda x: x.position);
	chained_events = [];
	for event in sorted_events:
		if (len(chained_events) > 0):
			last_size = chained_events[-1].size;
			last_endpoint = chained_events[-1].position + last_size;
			last_type = chained_events[-1].type;
			distance = abs(event.position - last_endpoint);
			if (distance <= (last_size + event.size) and last_type == event.type):
				chained_events[-1].close_cov = event.close_cov;
				chained_events[-1].max_cov = max(chained_events[-1].max_cov, event.max_cov);
				chained_events[-1].max_event_length = max(chained_events[-1].max_event_length, event.max_event_length);
				chained_events[-1].size = event.position + event.size - chained_events[-1].position + 1;
				# print 'Chaining: %s' % chained_events[-1].verbose_to_string();
			else:
				chained_events.append(event);
		else:
			chained_events.append(event);

	# sys.stderr.write('Chained events:\n');
	# for event in chained_events:
	# 	sys.stderr.write(event.csv_line() + '\n');

	return chained_events;
	# return events;



def calculate_stats(events, events_truth, fp_out):
	tp = 0;
	fp = 0;
	for event in events:
		if (event.evaluation == TP):
			tp += 1;
		elif (event.evaluation == FP):
			fp += 1;
	num_real_events = 0;
	for truth_event in events_truth:
		if ('real' in truth_event.comment):
			num_real_events += 1;

	recall = float(tp) / float(num_real_events);
	if (tp + fp) > 0:
		precision = float(tp) / float(tp + fp)
	else:
		precision = 0.0
	F1 = 0.0 if (tp == 0) else (2.0 * (precision * recall) / (precision + recall));

	fp_out.write('Precision: %.1f%% (%d / %d)\n' % (precision * 100.0, tp, (tp + fp)));
	fp_out.write('Recall: %.1f%% (%d / %d)\n' % (recall * 100.0, tp, (num_real_events)));
	fp_out.write('F1: %.1f%%\n' % (F1 * 100.0));



def main():

	if (len(sys.argv) < 3):
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s default\n' % (sys.argv[0]));
		sys.stderr.write('\tor\n');
		sys.stderr.write('\t%s <path_to_truths_csv> <events1.csv> [<events2.csv> ...]\n' % (sys.argv[0]));
		sys.stderr.write('\n');
		sys.stderr.write('Using default files for legacy.\n');
		sys.stderr.write('\n');
		exit(0);

	### This is the full list. Commented out when I tested on the latest version of the anchored alignment.
	# if (len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] == 'default')):
		# events_truth_path = 'results/truth_events-full.csv';
		# all_events_paths = ['results/events-indel_events-graphmap.csv',
		# 					'results/events-indel_events-bwamem.csv',
		# 					'results/events-indel_events-last.csv',
		# 					'results/events-indel_events-last-nanopore.csv',
		# 					'results/events-indel_events-blasr.csv'];
	# else:
	events_truth_path = sys.argv[1];
	all_events_paths = [];
	for event_file in sys.argv[2:]:
		all_events_paths.append(event_file);



	# all_events_paths = ['results/events-indel_events-graphmap-20150525.csv'];
###	all_events_paths = ['results/events-indel_events-graphmap-20150527.csv'];

	# all_events_paths = ['results/events-indel-graphmap.csv'];
	# all_events_paths = ['results/events-indel_events-blasr.csv'];

	# events_for_evaluation_path = 'results/events-indel-graphmap.csv';
	# events_for_evaluation_path = 'results/events-indel_events-bwamem.csv';
	# events_for_evaluation_path = 'results/events-indel_events-blasr.csv';
# znaci problem je u tome sto u ovim konkretnim primjerima postoji read koji ima event, ali jako veliki, preveliki od ocekivanoga
# vrlo vjerojatno lose filtrirano
# prije se to nije pojavljivalo
# max_event_length je preveliki
# pogledati figure za "channel_51_read_42_twodirections_/dev/shm/downloads/LomanLabz_E.coli_MG1655_3311_1_ch51_file42_strand.fast5"



	sys.stderr.write('Loading truth events...\n');
	events_truth = load_events(events_truth_path);

	for events_for_evaluation_path in all_events_paths:
		dirname_path = os.path.dirname(events_truth_path);
		if (len(dirname_path) == 0 or (dirname_path[0] in './') == False):
			dirname_path = './%s' % dirname_path;
		out_path = '%s/evaluated-%s' % (dirname_path, os.path.basename(events_for_evaluation_path));
		sys.stderr.write('Processing file "%s" and writing output to "%s"...\n' % (events_for_evaluation_path, out_path));
		fp_out = open(out_path, 'w');
		events = load_events(events_for_evaluation_path);
		evaluate(events, events_truth, fp_out);
		calculate_stats(events, events_truth, fp_out);
		sys.stderr.write('Stats for "%s":\n' % (os.path.basename(events_for_evaluation_path)));
		calculate_stats(events, events_truth, sys.stderr);
		sys.stderr.write('\n');
		fp_out.close();
	sys.stderr.write('Done!\n\n');

	# events_for_evaluation_path = 'results/events-blacklist-bwamem.csv';
	# events_for_evaluation_path = 'results/events-blacklist_events-blasr.csv';

# ./evaluate_events.py > results/test-BWA-MEM.csv
# ./evaluate_events.py > results/test-GraphMap.csv
# ./evaluate_events.py > results/test-BLASR.csv


	# for event in events:
	# 	event.verbose();
	# 	sys.stderr.write('\n');



if __name__ == "__main__":
	TP = 'true';
	FP = 'false';
	TN = 'blacklisted';

	main();

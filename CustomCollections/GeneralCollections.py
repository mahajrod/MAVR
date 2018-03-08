#!/usr/bin/env python

from collections import OrderedDict, MutableSet, Iterable


class TwoLvlDict(OrderedDict):

    def __init__(self, dictionary=None, input_file=None, column_sep="\t", value_sep=",",
                 split_values=False, absent_symbol=None, value_expression=None, ignore_value_repeats=False):
        OrderedDict.__init__(self)
        if dictionary:
            for fl_key in dictionary:
                self[fl_key] = OrderedDict()
                for sl_key in dictionary[fl_key]:
                    self[fl_key][sl_key] = dictionary[fl_key][sl_key]
        if input_file:
            input_file_list = [input_file] if isinstance(input_file, str) else input_file
            for filename in input_file_list:
                with open(filename, "r") as in_fd:
                    header_list = in_fd.readline().strip().split(column_sep)
                    for fl_key in header_list[1:]:
                        if fl_key not in self:
                            self[fl_key] = OrderedDict()
                    for line in in_fd:
                        tmp_list = line.strip().split("\t")
                        if len(tmp_list) != len(header_list):
                            raise ValueError("Error!!! Several values are absent! Malformed input file %s" % filename)
                        for i in range(1, len(header_list)): #in tmp_list[1:]:
                            #if tmp_list[0] not in self[header_list[i]]:
                            #    self[header_list[i]][tmp_list[0]] =
                            if not ignore_value_repeats:
                                if tmp_list[0] in self[header_list[i]]:
                                    if self[header_list[i]][tmp_list[0]] is not None:
                                        raise ValueError("Value is already present in TwoLvlDict. Filename: %s. Fl_key: %s. Sl_key: %s" % (filename, header_list[i], tmp_list[0]))

                            value = tmp_list[i]
                            if absent_symbol:
                                if value == absent_symbol:
                                    continue
                            if split_values:
                                value = value.split(value_sep)
                                if value_expression:
                                    value = map(value_expression, value)
                            else:
                                if value_expression:
                                    value = value_expression(value)
                            self[header_list[i]][tmp_list[0]] = value

    def table_form(self, absent_symbol="0", sort=True, column_sep="\t", list_sep=","):
        first_level_keys = list(self.keys())
        second_level_keys = OrderedSet()
        for fl_key in first_level_keys:
            for sl_key in self[fl_key]:
                second_level_keys.add(sl_key)
        if sort:
            first_level_keys.sort()
            second_level_keys = sorted(list(second_level_keys))
        string = "." + column_sep + column_sep.join([str(x) for x in first_level_keys]) + "\n"

        for sl_key in second_level_keys:
            key_counts_list = []
            for fl_key in first_level_keys:
                if sl_key not in self[fl_key]:
                    key_counts_list.append(absent_symbol)
                else:
                    list_or_tuple_type = isinstance(self[fl_key][sl_key], list) \
                                         or isinstance(self[fl_key][sl_key], tuple)
                    key_counts_list.append(list_sep.join(map(lambda y: str(y), self[fl_key][sl_key]))
                                           if list_or_tuple_type else str(self[fl_key][sl_key]))
            string += str(sl_key) + column_sep + column_sep.join(key_counts_list) + "\n"
        return string

    def sl_keys(self):
        sl_key_set = OrderedSet()
        for fl_key in self:
            for sl_key in self[fl_key]:
                sl_key_set.add(sl_key)
        return sl_key_set

    def all_values(self):
        # iterate over table by columns
        for fl_key in self:
            for sl_key in self[fl_key]:
                yield self[fl_key][sl_key]

    def filter_by_value(self, expression):
        for fl_key in self:
            for sl_key in self[fl_key]:
                if not expression(self[fl_key][sl_key]):
                    self[fl_key].pop(sl_key, None)
                    if not self[fl_key]:
                        self.pop(fl_key, None)

    def filter_by_line(self, expression):
        filtered_out = TwoLvlDict()
        for sl_key in self.sl_keys():
            line_list = []
            for fl_key in self:
                if sl_key in self[fl_key]:
                    line_list.append(self[fl_key][sl_key])
            if not expression(line_list):
                for fl_key in self:
                    if sl_key not in self[fl_key]:
                        continue
                    if fl_key not in filtered_out:
                        filtered_out[fl_key] = OrderedDict()
                    filtered_out[fl_key][sl_key] = self[fl_key].pop(sl_key, None)
        return filtered_out

    def write(self, out_filename, absent_symbol="0", close_after_if_file_object=False, sort=True):
        if isinstance(out_filename, file):
            out_filename.write(self.table_form(absent_symbol=absent_symbol, sort=sort))
            if close_after_if_file_object:
                out_filename.close()
        else:
            with open(out_filename, "w") as out_fd:
                out_fd.write(self.table_form(absent_symbol=absent_symbol, sort=sort))

    def write_splited(self, out_dir="./", extension="t", value_separator=","):
        from Routines.File import check_path
        for fl_key in self:
            with open("%s%s.%s" % (check_path(out_dir), fl_key, extension), "w") as out_fd:
                for sl_key in self[fl_key]:
                    out_fd.write("%s\t%s\n" % (sl_key, value_separator.join(self[fl_key][sl_key])))

    def read(self, filename, absent_symbol="0", column_sep="\t", value_handler=None):
        with open(filename, "r") as in_fd:
            fl_keys = in_fd.readline().strip().split(column_sep)
            for fl_key in fl_keys[1:]:
                self[fl_key] = OrderedDict()
            for line in in_fd:
                line_list = line.strip().split(column_sep)
                sl_key = line_list[0]
                for index in range(1, len(line_list)):
                    if line_list[index] != absent_symbol:
                        fl_key = fl_keys[index]
                        self[fl_key][sl_key] = value_handler[line_list[index]] if value_handler is not None else line_list[index]


class OrderedSet(MutableSet):

    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

    @staticmethod
    def union(*sets):
        union = OrderedSet()
        union.union(*sets)
        return union

    def union(self, *sets):
        for set in sets:
            self |= set


class IdList(list):

    def __init__(self, idlist=None, filename=None, header=False, close_after_if_file_object=False, column_number=None,
                 column_separator="\t", comments_prefix=None, id_in_column_separator=None):
        list.__init__(self)

        if filename:
            self.read(filename, header=header, close_after_if_file_object=close_after_if_file_object,
                      column_number=column_number, column_separator=column_separator,
                      comments_prefix=comments_prefix, id_in_column_separator=id_in_column_separator)
        elif idlist:
            for element in idlist:
                self.append(element)

    def read(self, filename, header=False, close_after_if_file_object=False, column_number=None, column_separator="\t",
             comments_prefix=None, id_in_column_separator=None):
        #reads ids from file with one id per line
        
        in_fd = filename if isinstance(filename, file) else open(filename, "r")
        if comments_prefix:
            com_pref_len = len(comments_prefix)
        if header:
            self.header = in_fd.readline().strip()
        for line in in_fd:
            if comments_prefix:
                if line[: com_pref_len] == comments_prefix:
                    continue

            column = line.strip().split(column_separator)[column_number] if column_number is not None else line.strip()
            if id_in_column_separator:
                self.extend(column.split(id_in_column_separator))
            else:
                self.append(column)
        if (not isinstance(filename, file)) or close_after_if_file_object:
            in_fd.close()
        return self
    
    def write(self, filename, header=False, close_after_if_file_object=False):
        out_fd = filename if isinstance(filename, file) else open(filename, "w")

        if header:
            if header is True and self.header:
                out_fd.write(self.header + "\n")
        for entry in self:
            out_fd.write(str(entry) + "\n")

        if (not isinstance(filename, file)) or close_after_if_file_object:
            out_fd.close()
            

class IdSet(OrderedSet):

    def __init__(self, idset=None, filename=None, header=False, close_after_if_file_object=False, column_number=None,
                 column_separator="\t", comments_prefix=None, id_in_column_separator=None):
        OrderedSet.__init__(self)

        if filename:
            self.read(filename, header=header, close_after_if_file_object=close_after_if_file_object,
                      column_number=column_number, column_separator=column_separator,
                      comments_prefix=comments_prefix, id_in_column_separator=id_in_column_separator)
        if idset:
            for element in idset:
                self.add(element)

    def read(self, filename, header=False, close_after_if_file_object=False, column_number=None, column_separator="\t",
             comments_prefix=None, id_in_column_separator=None):
        #reads ids from file with one id per line
        
        in_fd = filename if isinstance(filename, file) else open(filename, "r")
        if comments_prefix:
            com_pref_len = len(comments_prefix)
        if header:
            self.header = in_fd.readline().strip()
        for line in in_fd:
            if comments_prefix:
                if line[: com_pref_len] == comments_prefix:
                    continue
            ids = line.strip().split(column_separator)[column_number] if column_number is not None else line.strip()
            if id_in_column_separator:
                ids = set(ids.split(id_in_column_separator))
                for entry in ids:
                    self.add(entry)
            else:
                self.add(ids)

        if (not isinstance(filename, file)) or close_after_if_file_object:
            in_fd.close()
        return self

    def write(self, filename, header=False, close_after_if_file_object=False, sort=True):
        
        out_fd = filename if isinstance(filename, file) else open(filename, "w")
        if header:
            if header is True and self.header:
                out_fd.write(self.header + "\n")
        for entry in sorted(self) if sort else self:
            out_fd.write(entry + "\n")

        if (not isinstance(filename, file)) or close_after_if_file_object:
            out_fd.close()
            

class WDict(OrderedDict):
    def write(self, outfile, header=None, separator="\t", close_after_if_file_object=False):
        
        out_fd = outfile if isinstance(outfile, file) else open(outfile, "w")
        if header:
            if header is True and self.header:
                out_fd.write(self.header + "\n")
        for key in self:
            out_fd.write("%s%s%s\n" % (key, separator, self[key]))

        if (not isinstance(outfile, file)) or close_after_if_file_object:
            out_fd.close()


class SynDict(OrderedDict):

    def __init__(self, input=None, filename=None, header=False, separator="\t", allow_repeats_of_key=False,
                 split_values=False, values_separator=",", key_index=0, value_index=1,
                 close_after_if_file_object=False, expression=None, comments_prefix=None,
                 include_line_expression=None, include_value_expression=None):
        OrderedDict.__init__(self)

        if filename:
            self.read(filename, header=header, separator=separator,
                      allow_repeats_of_key=allow_repeats_of_key,
                      split_values=split_values, values_separator=values_separator,
                      key_index=key_index, value_index=value_index,
                      close_after_if_file_object=close_after_if_file_object, expression=expression,
                      comments_prefix=comments_prefix,
                      include_line_expression=include_line_expression,
                      include_value_expression=include_value_expression)
        elif input:
            OrderedDict.__init__(self, input)

    def count_synonyms(self):
        count_dict = OrderedDict()
        for key in self:
            if isinstance(self[key], Iterable) and (not isinstance(self[key], str)):
                count_dict[key] = len(self[key])
            else:
                count_dict[key] = 1

        return count_dict

    def read(self, filename, header=False, separator="\t", allow_repeats_of_key=False,
             split_values=False, values_separator=",", key_index=0, value_index=1,
             close_after_if_file_object=False, expression=None, comments_prefix=None,
             include_line_expression=None, include_value_expression=None):
        #reads ids from file with one id per line

        """
        IMPORTANT!!! Option allow_repeats_of_keys forces split_values.
        """
        # reads synonyms from file

        in_fd = filename if isinstance(filename, file) else open(filename, "r")
        if comments_prefix:
            com_pref_len = len(comments_prefix)

        if header:
            for line in in_fd:
                if comments_prefix:
                    if line[: com_pref_len] == comments_prefix:
                        continue

                    else:
                        self.header = line.strip()
                        break
        else:
            self.header = None

        for line in in_fd:
            if comments_prefix:
                if line[: com_pref_len] == comments_prefix:
                    continue
            if "\t0\t" in line:
                print line

            if include_line_expression:
                if not include_line_expression(line):
                    if "\t0\t" in line:
                        print line
                    continue
            #key, value = line.strip().split(separator)

            tmp = line.strip().split(separator) if separator else line.strip().split()
            #print tmp
            #print line
            #print tmp
            key, value = tmp[key_index], tmp[value_index] if isinstance(value_index, int) else [tmp[index_entry] for index_entry in value_index]
            #print key, value
            if split_values or allow_repeats_of_key:
                if isinstance(value, str):
                    value = value.split(values_separator)
                else:
                    tmp = []
                    for entry in value:
                        tmp += entry.split(values_separator)
                    value = tmp
            #print key, value
            if expression:
                if split_values or allow_repeats_of_key:
                    value = map(expression, value)
                else:
                    value = expression(value)

            if include_value_expression:
                if split_values or allow_repeats_of_key:
                    filtered_value_list = []
                    for entry in value:
                        if include_value_expression(entry):
                            filtered_value_list.append(entry)
                    if filtered_value_list:
                        value = filtered_value_list
                    else:
                        continue
                else:
                    if not include_value_expression(value):
                        print value
                        continue

            if key not in self:
                self[key] = value
            else:
                if allow_repeats_of_key:
                    self[key] += value
                else:
                    print "Repeated key: %s" % key
                    print self[key]
                    raise ValueError("Error while reading to SynDit: key is repeated")

        if (not isinstance(filename, file)) or close_after_if_file_object:
            in_fd.close()
        return self

    def count_all_synonyms(self):
        number = 0

        for key in self:
            if isinstance(self[key], Iterable) and (not isinstance(self[key], str)):
                number += len(self[key])
            else:
                number += 1

        return number

    def write(self, filename, header=False, separator="\t",
              splited_values=False, values_separator=",",
              close_after_if_file_object=False):

        # reads synonyms from file
        out_fd = filename if isinstance(filename, file) else open(filename, "w")
        if header:
            if header:
                if header is True and self.header:
                    out_fd.write(self.header + "\n")

        for entry in self:
            #print "AAAAAAAAAAAAAA"
            out_fd.write("%s%s%s\n" % (entry, separator,
                                       values_separator.join(map(str, self[entry])) if splited_values else str(self[entry])))
        if (not isinstance(filename, file)) or close_after_if_file_object:
            out_fd.close()

    def remove_value_repeats(self):
        collapsed_dict = SynDict()
        for key in self:
            collapsed_dict[key] = set(self[key])
        return collapsed_dict

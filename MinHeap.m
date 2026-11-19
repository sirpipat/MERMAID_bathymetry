classdef MinHeap < handle
%MINHEAP  General-purpose binary min-heap (priority queue).
%
%   Each entry is a pair (key, item):
%     key   numeric row vector of fixed length K. Keys are compared
%           lexicographically: first element first; only if equal, the
%           second; and so on. Fully equal keys come out in insertion order.
%     item  anything at all (struct, cell, number, object). The heap never
%           looks inside it, so new fields can be added without touching
%           this class.
%
%   Usage
%     h = MinHeap(2);                       % keys have 2 elements
%     h.push([5 3], struct('name','x'));
%     h.push([5 1], struct('name','y'));
%     [key, item] = h.pop();                % key = [5 1], item.name = 'y'
%     h.isEmpty()                           % false
%     h.Count                               % 1
%
%   MinHeap is a handle class: push/pop change the heap in place, so there
%   is no need to write h = h.push(...).
%
% Last modified by spipatprathanporn@ucsd.edu, 11/18/2026

    properties (SetAccess = private)
        Count = 0       % number of entries currently stored
    end

    properties (Access = private)
        Keys            % capacity-by-K matrix, row i = key of slot i
        Seq             % capacity-by-1, insertion counter for tie-breaking
        Items           % capacity-by-1 cell, stored items
        NextSeq = 1     % next insertion counter value
    end

    methods
        function obj = MinHeap(keyLength, capacity)
        %MINHEAP  Construct an empty heap.
        %
        %   INPUTS
        %     keyLength  number of elements in every key (K >= 1)
        %     capacity   (optional) initial number of slots, default 64.
        %                Storage doubles automatically when full.
        %
        %   OUTPUTS
        %     obj        the new, empty MinHeap

            if nargin < 2 || isempty(capacity)
                capacity = 64;
            end
            obj.Keys  = zeros(capacity, keyLength);
            obj.Seq   = zeros(capacity, 1);
            obj.Items = cell(capacity, 1);
        end

        function push(obj, key, item)
        %PUSH  Insert an item with the given key.
        %
        %   INPUTS
        %     obj    the heap
        %     key    numeric vector with keyLength elements; smaller
        %            (lexicographically) is popped first
        %     item   any value to store with the key
        %
        %   OUTPUTS
        %     none (the heap is modified in place)

            key = key(:).';
            if numel(key) ~= size(obj.Keys, 2)
                error('MinHeap:keyLength', ...
                      'Key must have %d elements.', size(obj.Keys, 2));
            end
            if obj.Count == size(obj.Keys, 1)       % grow storage
                cap = 2 * size(obj.Keys, 1);
                obj.Keys(cap, :)  = 0;
                obj.Seq(cap, 1)   = 0;
                obj.Items{cap, 1} = [];
            end

            obj.Count = obj.Count + 1;
            i = obj.Count;
            obj.Keys(i, :) = key;
            obj.Seq(i)     = obj.NextSeq;
            obj.Items{i}   = item;
            obj.NextSeq    = obj.NextSeq + 1;
            obj.siftUp(i);
        end

        function [key, item] = pop(obj)
        %POP  Remove and return the entry with the smallest key.
        %
        %   INPUTS
        %     obj    the heap (must not be empty)
        %
        %   OUTPUTS
        %     key    key of the removed entry (row vector)
        %     item   the stored item

            if obj.Count == 0
                error('MinHeap:empty', 'Cannot pop from an empty heap.');
            end
            key  = obj.Keys(1, :);
            item = obj.Items{1};

            last = obj.Count;
            obj.Keys(1, :) = obj.Keys(last, :);
            obj.Seq(1)     = obj.Seq(last);
            obj.Items{1}   = obj.Items{last};
            obj.Items{last} = [];                   % release memory
            obj.Count = last - 1;
            obj.siftDown(1);
        end

        function [key, item] = peek(obj)
        %PEEK  Return the smallest entry without removing it.
        %
        %   INPUTS
        %     obj    the heap (must not be empty)
        %
        %   OUTPUTS
        %     key    key of the smallest entry
        %     item   its item

            if obj.Count == 0
                error('MinHeap:empty', 'Cannot peek into an empty heap.');
            end
            key  = obj.Keys(1, :);
            item = obj.Items{1};
        end

        function tf = isEmpty(obj)
        %ISEMPTY  True if the heap has no entries.
        %
        %   INPUTS
        %     obj    the heap
        %
        %   OUTPUTS
        %     tf     logical scalar

            tf = (obj.Count == 0);
        end
    end

    methods (Access = private)
        function siftUp(obj, i)
        %SIFTUP  Move slot i up until its parent is not larger.
        %
        %   INPUTS
        %     obj    the heap
        %     i      slot index of a newly inserted entry
        %
        %   OUTPUTS
        %     none (the heap is modified in place)

            while i > 1
                parent = floor(i / 2);
                if obj.less(i, parent)
                    obj.swap(i, parent);
                    i = parent;
                else
                    break;
                end
            end
        end

        function siftDown(obj, i)
        %SIFTDOWN  Move slot i down until both children are not smaller.
        %
        %   INPUTS
        %     obj    the heap
        %     i      slot index (1 after a pop)
        %
        %   OUTPUTS
        %     none (the heap is modified in place)

            n = obj.Count;
            while true
                l = 2 * i;
                r = l + 1;
                s = i;
                if l <= n && obj.less(l, s), s = l; end
                if r <= n && obj.less(r, s), s = r; end
                if s == i, break; end
                obj.swap(i, s);
                i = s;
            end
        end

        function tf = less(obj, i, j)
        %LESS  True if slot i should be popped before slot j.
        %
        %   INPUTS
        %     obj    the heap
        %     i, j   slot indices
        %
        %   OUTPUTS
        %     tf     logical: lexicographic key comparison, then insertion order

            ki = obj.Keys(i, :);
            kj = obj.Keys(j, :);
            d = find(ki ~= kj, 1);
            if isempty(d)
                tf = obj.Seq(i) < obj.Seq(j);
            else
                tf = ki(d) < kj(d);
            end
        end

        function swap(obj, i, j)
        %SWAP  Exchange the contents of slots i and j.
        %
        %   INPUTS
        %     obj    the heap
        %     i, j   slot indices
        %
        %   OUTPUTS
        %     none (the heap is modified in place)

            obj.Keys([i j], :) = obj.Keys([j i], :);
            obj.Seq([i j])     = obj.Seq([j i]);
            obj.Items([i j])   = obj.Items([j i]);
        end
    end
end
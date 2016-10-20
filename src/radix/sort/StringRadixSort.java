/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package radix.sort;

/**
 * Radix sort for byte strings
 */
public class StringRadixSort {

	private static final int INSERTION_SORT_THRESHOLD = 64;
	private static final int BITS_PER_LEVEL = 8;
	private static final int NODE_SIZE_CHAR = 1 << BITS_PER_LEVEL;

	private StringRadixSort() {
		throw new AssertionError();
	}

	/**
	 * @param data    strings, serialized as byte, zero terminated
	 * @param offsets
	 * @param size
	 * @return
	 */
	public static int[] radixSortWithPermutationTrace(byte[] data, int[] offsets, int size) {
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortStringsWithPermutationVector(data, offsets, new long[offsets.length], 0, size, 0, permutationVector);
		return permutationVector;
	}

	/**
	 * Sort strings (bytes, zero terminated) indirectly by sorting their staring offsets.
	 * Permutation vector keeps track about how the order was changed with respect to the original order.
	 * <p>
	 * TODO: if you do not need the changes with respect to the original data, remove the permutation vector stuff...
	 *
	 * @param data              byte array of multiple zero terminated strings (as bytes as in ASCII or UTF-8
	 * @param offsets           array with the start offsets of all strings in data
	 * @param keys              helper structure for efficiency
	 * @param startOffset
	 * @param end
	 * @param currentChar
	 * @param permutationVector here we keep track how the strings have been reordered with respect to original order.
	 */
	private static void radixSortStringsWithPermutationVector(final byte[] data, int[] offsets, final long[] keys, int startOffset, final int end, final int currentChar, int[] permutationVector) {
		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE_CHAR];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE_CHAR];
		//emulates currentChar % 8
		final int longOffset = currentChar & 7;
		final int currentShift = (7 - longOffset) << 3;
		//load needed keys
		if (longOffset == 0) {
			boolean canSkip = startOffset < end;
			for (int i = startOffset; i < end; ++i) {
				long key = BitUtils.extractBytesIntoLong(data, offsets[i] + currentChar);
				keys[i] = key;
				canSkip &= (keys[startOffset] == key);
			}

			if (canSkip) {
				if (!BitUtils.containsZeroByte(keys[startOffset])) {//we do not recurs if the string was terminated
					radixSortStringsWithPermutationVector(data, offsets, keys, startOffset, end, currentChar + 8, permutationVector);
				}
				return;
			}
		}

		//obtain bucket counts
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[(int) ((keys[i] >> currentShift) & 0xFF)];
		}


		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < bucketCounts.length; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		for (int i = 0; i < bucketCounts.length; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (i > 0) {//do not recurs for the first bucket, because it contains only ended strings
						radixSortStringsWithPermutationVector(data, offsets, keys, startOffset, end, currentChar + 1, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//for all buckets:
		for (int i = NODE_SIZE_CHAR; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curItem = offsets[bucketPointers[i]];
				long curKey = keys[bucketPointers[i]];
				int curIdx = permutationVector[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = (int) ((curKey >> currentShift) & 0xFF);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpItem = offsets[bucketPointers[itemKey]];
					long tmpKey = keys[bucketPointers[itemKey]];
					int tmpIdx = permutationVector[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					offsets[bucketPointers[itemKey]] = curItem;
					keys[bucketPointers[itemKey]] = curKey;
					permutationVector[bucketPointers[itemKey]++] = curIdx;

					//make the removed item to curItem
					curItem = tmpItem;
					curKey = tmpKey;
					curIdx = tmpIdx;
					//calculate the key for this item and continue until we find an item belonging to bucket i

					itemKey = (int) ((curKey >> currentShift) & 0xFF);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				offsets[bucketPointers[i]] = curItem;
				keys[bucketPointers[i]] = curKey;
				permutationVector[bucketPointers[i]++] = curIdx;

			}
		}

		//for all buckets (except the very first, containing only ended strings)...
		for (int i = NODE_SIZE_CHAR; --i >= 1; ) {
			final int currentPointer = bucketPointers[i];
			//calculate the bucket size
			final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
			//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
			if (bucketSize > INSERTION_SORT_THRESHOLD) {
				radixSortStringsWithPermutationVector(data, offsets, keys, currentPointer - bucketSize, currentPointer, currentChar + 1, permutationVector);
			} else if (bucketSize > 1) {
				insertionSortWithPermutationVector(data, offsets, keys, currentPointer - bucketSize, currentPointer, currentChar + 1, permutationVector);
			}
		}
	}

	private static int findSkipPotential(long and, long or) {
		return Long.numberOfLeadingZeros(and ^ or) >> 3;
	}

	/**
	 * @param data
	 * @param offset
	 * @param end
	 */
	private static void insertionSortWithPermutationVector(byte[] data, int[] offsets, long[] keys, int offset, int end, int curByte, int[] permutationVector) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && (keys[j - 1] > keys[j] || keys[j - 1] == keys[j] && (compare(data, offsets[j - 1] + curByte, offsets[j] + curByte) > 0)); --j) {
				final int tmpItem = offsets[j];
				final long tmpKey = keys[j];
				final int tmpIdx = permutationVector[j];
				keys[j] = keys[j - 1];
				keys[j - 1] = tmpKey;
				offsets[j] = offsets[j - 1];
				offsets[j - 1] = tmpItem;
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = tmpIdx;
			}
		}
	}

	private static int compare(final byte[] data, int off1, int off2) {
		while (data[off1] != 0) {
			int dif = (data[off1++] & 0xFF) - (data[off2++] & 0xFF);
			if (dif != 0) {
				return dif;
			}

		}

		return data[off1] - data[off2];
	}

	/**
	 * Generate fresh permutation vector to be used in the radix sorts that
	 * produce permutation traces.
	 *
	 * @param size
	 * @return int array of given size with array[i] == i.
	 */
	private static int[] generateFreshPermutationVector(int size) {
		final int[] permutationVector = new int[size];
		for (int i = 0; i < permutationVector.length; ++i) {
			permutationVector[i] = i;
		}
		return permutationVector;
	}

}

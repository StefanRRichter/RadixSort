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
 * Tuned implementation of Radix Sort for primitives.
 * <p>
 * This code has some crazy duplication for performance :-(
 */
public final class RadixSort {

	private static final int INSERTION_SORT_THRESHOLD = 64;
	private static final int BITS_PER_LEVEL = 8;
	private static final int NODE_SIZE = 1 << BITS_PER_LEVEL;
	private static final int NODE_SIZE_CHAR = 1 << BITS_PER_LEVEL;
	private static final int MASK = (1 << BITS_PER_LEVEL) - 1;

	private RadixSort() {
		throw new AssertionError();
	}

//-------------------------------PUBLIC INTERFACE-------------------------------

	/**
	 * Radix Sort for ints.
	 *
	 * @param array
	 */
	public static void radixSort(int[] array) {
		radixSort(array, 0, array.length, Integer.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for longs.
	 *
	 * @param array
	 */
	public static void radixSort(long[] array) {
		radixSort(array, 0, array.length, Long.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for floats.
	 *
	 * @param array
	 */
	public static void radixSort(float[] array) {
		radixSort(array, 0, array.length, Float.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for doubles.
	 *
	 * @param array
	 */
	public static void radixSort(double[] array) {
		radixSort(array, 0, array.length, Double.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for doubles.
	 *
	 * @param array
	 */
	public static void radixSort(byte[][] array) {
		radixSort(array, array.length);
	}

	/**
	 * Radix Sort for ints.
	 *
	 * @param array
	 */

	public static void radixSort(int[] array, int size) {
		radixSort(array, 0, size, Integer.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for longs.
	 *
	 * @param array
	 */
	public static void radixSort(long[] array, int size) {
		radixSort(array, 0, size, Long.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for floats.
	 *
	 * @param array
	 */
	public static void radixSort(float[] array, int size) {
		radixSort(array, 0, size, Float.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for doubles.
	 *
	 * @param array
	 */
	public static void radixSort(double[] array, int size) {
		radixSort(array, 0, size, Double.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for doubles.
	 *
	 * @param array
	 */
	public static void radixSort(byte[][] array, int size) {
		radixSort(array, 0, size, 0);
	}

//------------------------------------------------------------------------------

	/**
	 * Radix Sort for unsigned ints.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(int[] array) {
		radixSortUnsigned(array, 0, array.length, Integer.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned longs.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(long[] array) {
		radixSortUnsigned(array, 0, array.length, Long.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned floats.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(float[] array) {
		radixSortUnsigned(array, 0, array.length, Float.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned doubles.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(double[] array) {
		radixSortUnsigned(array, 0, array.length, Double.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned ints.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(int[] array, int size) {
		radixSortUnsigned(array, 0, size, Integer.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned longs.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(long[] array, int size) {
		radixSortUnsigned(array, 0, size, Long.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned floats.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(float[] array, int size) {
		radixSortUnsigned(array, 0, size, Float.SIZE - BITS_PER_LEVEL);
	}

	/**
	 * Radix Sort for unsigned doubles.
	 * <p/>
	 * Faster but does not handle negative values correctly.
	 *
	 * @param array
	 */
	public static void radixSortUnsigned(double[] array, int size) {
		radixSortUnsigned(array, 0, size, Double.SIZE - BITS_PER_LEVEL);
	}

//------------------------------------------------------------------------------

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(int[] array) {
		return radixSortWithPermutationTrace(array, array.length);
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(long[] array) {
		return radixSortWithPermutationTrace(array, array.length);
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(float[] array) {
		return radixSortWithPermutationTrace(array, array.length);
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(double[] array) {
		return radixSortWithPermutationTrace(array, array.length);
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(byte[][] array) {
		return radixSortWithPermutationTrace(array, array.length);
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(int[] array, int size) {
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortWithPermutationTrace(array, 0, size, Integer.SIZE - BITS_PER_LEVEL, permutationVector);
		return permutationVector;
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(long[] array, int size) {
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortWithPermutationTrace(array, 0, size, Long.SIZE - BITS_PER_LEVEL, permutationVector);
		return permutationVector;
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(float[] array, int size) {
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortWithPermutationTrace(array, 0, size, Float.SIZE - BITS_PER_LEVEL, permutationVector);
		return permutationVector;
	}

	/**
	 * Radix Sort that keeps track on the element movement (so that you could
	 * reorder other arrays accordingly). Of course this has some performance
	 * penalty.
	 *
	 * @param array
	 * @return permutaton vector (from-to mapping showing how elements were moved)
	 */
	public static int[] radixSortWithPermutationTrace(double[] array, int size) {
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortWithPermutationTrace(array, 0, size, Double.SIZE - BITS_PER_LEVEL, permutationVector);
		return permutationVector;
	}

	/**
	 * @param array
	 * @param size
	 * @return
	 */
	public static int[] radixSortWithPermutationTrace(byte[][] array, int size) {
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortStringsWithPermutationVector(array, new long[array.length], 0, size, 0, permutationVector);
		return permutationVector;
	}

//---------------------------------RADIX SORTS----------------------------------

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSort(final int[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			//bucketkey = invert first bit (to handle negative values), shift for radix step and apply mask
			++bucketCounts[MASK & ((array[i] ^ (1 << Integer.SIZE - 1)) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSort(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & ((curItem ^ (1 << Integer.SIZE - 1)) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & ((curItem ^ (1 << Integer.SIZE - 1)) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSort(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSort(final long[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			//bucketkey = invert first bit (to handle negative values), shift for radix step and apply mask
			++bucketCounts[MASK & (int) ((array[i] ^ (1L << Long.SIZE - 1)) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSort(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				long curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (int) ((curItem ^ (1L << Long.SIZE - 1)) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					long tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (int) ((curItem ^ (1L << Long.SIZE - 1)) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSort(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSort(final float[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			//needed in java :-(
			final int itemKey = Float.floatToRawIntBits(array[i]);
			//bucketkey = invert first bit and invert if first bit was set (to handle negative floating point values)
			//then shift for radix step and apply mask
			++bucketCounts[MASK & ((itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1))) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSort(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				float curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = Float.floatToRawIntBits(curItem);
				itemKey = MASK & ((itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1))) >> bitShift);
//                int itemKey = curItem < 0f ? MASK & ((-Float.floatToRawIntBits(curItem)) >> bitShift) : MASK & ((Float.floatToRawIntBits(curItem) ^ (1 << Float.SIZE - 1)) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					float tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = Float.floatToRawIntBits(curItem);
					itemKey = MASK & ((itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1))) >> bitShift);
//                    itemKey = curItem < 0f ? MASK & ((-Float.floatToRawIntBits(curItem)) >> bitShift) : MASK & ((Float.floatToRawIntBits(curItem) ^ (1 << Float.SIZE - 1)) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSort(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSort(final double[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			//needed in java :-(
			final long itemKey = Double.doubleToRawLongBits(array[i]);
			//bucketkey = invert first bit and invert if first bit was set (to handle negative floating point values)
			//then shift for radix step and apply mask
			++bucketCounts[MASK & (int) ((itemKey ^ ((itemKey >> Double.SIZE - 1) | (1L << Double.SIZE - 1))) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSort(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				double curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				long curItemLong = Double.doubleToRawLongBits(curItem);
				int itemKey = MASK & (int) ((curItemLong ^ ((curItemLong >> Double.SIZE - 1) | (1L << Double.SIZE - 1))) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					double tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					curItemLong = Double.doubleToRawLongBits(curItem);
					itemKey = MASK & (int) ((curItemLong ^ ((curItemLong >> Double.SIZE - 1) | (1L << Double.SIZE - 1))) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSort(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}


	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param currentByte
	 */
	private static void radixSort(final byte[][] array, int startOffset, final int end, final int currentByte) {
		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		for (int i = startOffset; i < end; ++i) {
			byte[] current = array[i];
			if (current.length > currentByte) {
				++bucketCounts[current[currentByte] & 0xFF];
			} else {
				array[i] = array[startOffset];
				array[startOffset] = current;
				++startOffset;
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					radixSort(array, startOffset, end, currentByte + 1);

					return;
				} else {
					break;
				}
			}
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				byte[] curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = curItem[currentByte] & 0xFF;
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					byte[] tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;

					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i

					itemKey = curItem[currentByte] & 0xFF;
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;

			}
		}


		//if we are not in the final iteration...

		//calculate bitshift for the next iteration
		final int nextByte = currentByte + 1;
		//for all buckets...
		for (int i = NODE_SIZE; --i >= 0; ) {
			final int currentPointer = bucketPointers[i];
			//calculate the bucket size
			final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
			//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
			if (bucketSize > INSERTION_SORT_THRESHOLD) {
				radixSort(array, currentPointer - bucketSize, currentPointer, nextByte);
			} else if (bucketSize > 1) {
				insertionSort(array, currentPointer - bucketSize, currentPointer, currentByte);
			}
		}
	}

//-----------------------------UNSIGNED RADIX SORTS-----------------------------

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSortUnsigned(final int[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (array[i] >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortUnsigned(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (curItem >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (curItem >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortUnsigned(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSortUnsigned(final long[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (int) (array[i] >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortUnsigned(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				long curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (int) (curItem >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					long tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (int) (curItem >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortUnsigned(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSortUnsigned(final float[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (Float.floatToIntBits(array[i]) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortUnsigned(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				float curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (Float.floatToIntBits(curItem) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					float tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (Float.floatToIntBits(curItem) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > 1) {
					radixSortUnsigned(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 */
	private static void radixSortUnsigned(final double[] array, final int startOffset, final int end, final int bitShift) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (int) (Double.doubleToLongBits(array[i]) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortUnsigned(array, startOffset, end, bitShift - BITS_PER_LEVEL);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				double curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (int) (Double.doubleToLongBits(curItem) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					double tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (int) (Double.doubleToLongBits(curItem) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortUnsigned(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift);
				} else if (bucketSize > 1) {
					insertionSort(array, currentPointer - bucketSize, currentPointer);
				}
			}
		}
	}

//---------------------RADIX SORTS WITH PERMUTATION VECTOR----------------------

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 * @param permutationVector
	 */
	private static void radixSortWithPermutationTrace(final int[] array, final int startOffset, final int end, final int bitShift, final int[] permutationVector) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & ((array[i] ^ (1 << Integer.SIZE - 1)) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortWithPermutationTrace(array, startOffset, end, bitShift - BITS_PER_LEVEL, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curPerm = permutationVector[bucketPointers[i]];
				int curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & ((curItem ^ (1 << Integer.SIZE - 1)) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpPerm = permutationVector[bucketPointers[itemKey]];
					int tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curPerm;
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					curPerm = tmpPerm;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & ((curItem ^ (1 << Integer.SIZE - 1)) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curPerm;
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift, permutationVector);
				} else if (bucketSize > 1) {
					insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 * @param permutationVector
	 */
	private static void radixSortWithPermutationTrace(final long[] array, final int startOffset, final int end, final int bitShift, final int[] permutationVector) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (int) ((array[i] ^ (1L << Long.SIZE - 1)) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortWithPermutationTrace(array, startOffset, end, bitShift - BITS_PER_LEVEL, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curPerm = permutationVector[bucketPointers[i]];
				long curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (int) ((curItem ^ (1L << Long.SIZE - 1)) >> bitShift);

				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpPerm = permutationVector[bucketPointers[itemKey]];
					long tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curPerm;
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					curPerm = tmpPerm;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (int) ((curItem ^ (1L << Long.SIZE - 1)) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curPerm;
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift, permutationVector);
				} else if (bucketSize > 1) {
					insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 * @param permutationVector
	 */
	private static void radixSortWithPermutationTrace(final float[] array, final int startOffset, final int end, final int bitShift, final int[] permutationVector) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			final int itemKey = Float.floatToRawIntBits(array[i]);
			++bucketCounts[MASK & ((itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1))) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortWithPermutationTrace(array, startOffset, end, bitShift - BITS_PER_LEVEL, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curPerm = permutationVector[bucketPointers[i]];
				float curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = Float.floatToRawIntBits(curItem);
				itemKey = MASK & ((itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1))) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpPerm = permutationVector[bucketPointers[itemKey]];
					float tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curPerm;
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					curPerm = tmpPerm;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = Float.floatToRawIntBits(curItem);
					itemKey = MASK & ((itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1))) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curPerm;
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift, permutationVector);
				} else if (bucketSize > 1) {
					insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 * @param permutationVector
	 */
	private static void radixSortWithPermutationTrace(final double[] array, final int startOffset, final int end, final int bitShift, final int[] permutationVector) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			final long itemKey = Double.doubleToRawLongBits(array[i]);
			++bucketCounts[MASK & (int) ((itemKey ^ ((itemKey >> Double.SIZE - 1) | (1L << Double.SIZE - 1))) >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortWithPermutationTrace(array, startOffset, end, bitShift - BITS_PER_LEVEL, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curPerm = permutationVector[bucketPointers[i]];
				double curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				long curItemLong = Double.doubleToRawLongBits(curItem);
				int itemKey = MASK & (int) ((curItemLong ^ ((curItemLong >> Double.SIZE - 1) | (1L << Double.SIZE - 1))) >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpPerm = permutationVector[bucketPointers[itemKey]];
					double tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curPerm;
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					curPerm = tmpPerm;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					curItemLong = Double.doubleToRawLongBits(curItem);
					itemKey = MASK & (int) ((curItemLong ^ ((curItemLong >> Double.SIZE - 1) | (1L << Double.SIZE - 1))) >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curPerm;
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift, permutationVector);
				} else if (bucketSize > 1) {
					insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector);
				}
			}
		}
	}


	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param currentChar
	 */
	private static void radixSortStringsWithPermutationVector(final byte[][] array, final long[] keys, int startOffset, final int end, final int currentChar, int[] permutationVector) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE_CHAR];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE_CHAR];
		//emulates currentChar % 8
		final int longOffset = currentChar & 7;
		final int currentShift = (7 - longOffset) * 8;

		//Load needed keys
		if (longOffset == 0) {
			long and = 0xFFFFFFFFFFFFFFFFL;
			long or = 0;
			for (int i = end; --i >= startOffset; ) {
				long key = 0;
				final byte[] val = array[i];
				final int k = Math.min(val.length, currentChar + 8);
				for (int j = currentChar; j < k; ++j) // assume chars to be ascii
				{
					key = key << 8 | (long) val[j];
				}
				key = key << ((8 - k + currentChar) * 8);
				keys[i] = key;
				and &= key;
				or |= key;
			}
			if (and == or) {
				radixSortStringsWithPermutationVector(array, keys, startOffset, end, currentChar + 8, permutationVector);
				return;
			}
		}

		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[(int) ((keys[i] >> currentShift) & 0xFF)];
		}


		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (i > 0) {
						radixSortStringsWithPermutationVector(array, keys, startOffset, end, currentChar + 1, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//for all buckets:
		for (int i = NODE_SIZE_CHAR; --i >= 1; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				byte[] curItem = array[bucketPointers[i]];
				long curKey = keys[bucketPointers[i]];
				int curIdx = permutationVector[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = (int) ((curKey >> currentShift) & 0xFF);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					byte[] tmpItem = array[bucketPointers[itemKey]];
					long tmpKey = keys[bucketPointers[itemKey]];
					int tmpIdx = permutationVector[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					array[bucketPointers[itemKey]] = curItem;
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
				array[bucketPointers[i]] = curItem;
				keys[bucketPointers[i]] = curKey;
				permutationVector[bucketPointers[i]++] = curIdx;

			}
		}

		//for all buckets...
		for (int i = NODE_SIZE_CHAR; --i >= 1; ) {
			final int currentPointer = bucketPointers[i];
			//calculate the bucket size
			final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
			//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
			if (bucketSize > INSERTION_SORT_THRESHOLD) {
				radixSortStringsWithPermutationVector(array, keys, currentPointer - bucketSize, currentPointer, currentChar + 1, permutationVector);
			} else if (bucketSize > 1) {
				insertionSortWithPermutationVector(array, keys, currentPointer - bucketSize, currentPointer, currentChar + 1, permutationVector);
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param currentByte
	 * @param permutationVector
	 */
	@Deprecated
	private static void radixSortWithPermutationTrace(final byte[][] array, int startOffset, final int end, final int currentByte, final int[] permutationVector) {

		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		for (int i = startOffset; i < end; ++i) {
			byte[] current = array[i];
			if (current.length <= currentByte) {
				int tmp = permutationVector[i];
				array[i] = array[startOffset];
				permutationVector[i] = permutationVector[startOffset];
				array[startOffset] = current;
				permutationVector[startOffset] = tmp;
				++startOffset;
			} else {
				++bucketCounts[current[currentByte]];
			}
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (i > 0) {
						radixSortWithPermutationTrace(array, startOffset, end, currentByte + 1, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (i > 0) {
						radixSortWithPermutationTrace(array, startOffset, end, currentByte + 1, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				byte[] curItem = array[bucketPointers[i]];
				int curMapping = permutationVector[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = curItem[currentByte];
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					byte[] tmpItem = array[bucketPointers[itemKey]];
					int tmpMapping = permutationVector[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curMapping;
					array[bucketPointers[itemKey]++] = curItem;

					//make the removed item to curItem
					curItem = tmpItem;
					curMapping = tmpMapping;
					//calculate the key for this item and continue until we find an item belonging to bucket i

					itemKey = curItem[currentByte];
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curMapping;
				array[bucketPointers[i]++] = curItem;

			}
		}


		//if we are not in the final iteration...

		//calculate bitshift for the next iteration
		final int nextByte = currentByte + 1;
//            nextIterationBitShift = Math.max(0, nextIterationBitShift);
		//for all buckets...
		for (int i = NODE_SIZE; --i >= 1; ) {
			final int currentPointer = bucketPointers[i];
			//calculate the bucket size
			final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
			//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
			if (bucketSize > INSERTION_SORT_THRESHOLD) {
				radixSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextByte, permutationVector);
			} else if (bucketSize > 1) {
				insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector, currentByte);
			}
		}
	}

//--------------------------INSERTION SORTS FROM JDK----------------------------

	/**
	 * @param array
	 * @param offset
	 * @param end
	 */
	private static void insertionSort(int array[], int offset, int end) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				int temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 */
	private static void insertionSort(float array[], int offset, int end) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				float temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 */
	private static void insertionSort(long array[], int offset, int end) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				long temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 */
	private static void insertionSort(double array[], int offset, int end) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				double temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
			}
		}
	}


	/**
	 * @param array
	 * @param offset
	 * @param end
	 */
	private static void insertionSort(byte[][] array, int offset, int end, int curByte) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && (compare(array[j - 1], array[j], curByte) > 0); --j) {
				byte[] temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
			}
		}
	}

//-------------------INSERTION SORTS WITH PERMUTATION VECTOR--------------------

	/**
	 * @param array
	 * @param offset
	 * @param end
	 * @param permutationVector
	 */
	private static void insertionSortWithPermutationTrace(int array[], int offset, int end, int[] permutationVector) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				int temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
				int pmTemp = permutationVector[j];
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = pmTemp;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 * @param permutationVector
	 */
	private static void insertionSortWithPermutationTrace(long array[], int offset, int end, int[] permutationVector) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				long temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
				int pmTemp = permutationVector[j];
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = pmTemp;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 * @param permutationVector
	 */
	private static void insertionSortWithPermutationTrace(float array[], int offset, int end, int[] permutationVector) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				float temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
				int pmTemp = permutationVector[j];
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = pmTemp;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 * @param permutationVector
	 */
	private static void insertionSortWithPermutationTrace(double array[], int offset, int end, int[] permutationVector) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && array[j - 1] > array[j]; --j) {
				double temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
				int pmTemp = permutationVector[j];
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = pmTemp;
			}
		}
	}

	private static int compare(final byte[] bb1, final byte[] bb2, int curByte) {
		final int len = Math.min(bb1.length, bb2.length);
		while (curByte < len) {
			int dif = (bb1[curByte] & 0xFF) - (bb2[curByte] & 0xFF);
			if (dif != 0) {
				return dif;
			}
			++curByte;
		}
		return bb1.length - bb2.length;
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 */
	private static void insertionSortWithPermutationVector(byte[][] array, long[] keys, int offset, int end, int curByte, int[] permutationVector) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && (keys[j - 1] > keys[j] || keys[j - 1] == keys[j] && (compare(array[j - 1], array[j], curByte) > 0)); --j) {
				final byte[] tmpItem = array[j];
				final long tmpKey = keys[j];
				final int tmpIdx = permutationVector[j];
				keys[j] = keys[j - 1];
				keys[j - 1] = tmpKey;
				array[j] = array[j - 1];
				array[j - 1] = tmpItem;
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = tmpIdx;
			}
		}
	}

	/**
	 * @param array
	 * @param offset
	 * @param end
	 * @param permutationVector
	 */
	@Deprecated
	private static void insertionSortWithPermutationTrace(byte[][] array, int offset, int end, int[] permutationVector, int curByte) {
		for (int i = offset; i < end; ++i) {
			for (int j = i; j > offset && (compare(array[j - 1], array[j], curByte) > 0); --j) {
				byte[] temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
				int pmTemp = permutationVector[j];
				permutationVector[j] = permutationVector[j - 1];
				permutationVector[j - 1] = pmTemp;
			}
		}
	}

//-------------------GENERATION OF FRESH PERMUTATION VECTORS--------------------

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


	//-------------------------------


	public static void encode(int[] data, int start, int end) {
		for (int i = start; i < end; ++i) {
			final int itemKey = data[i];
			data[i] = (itemKey ^ ((itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1)));
		}
	}

	public static void decode(int[] data, int start, int end) {
		for (int i = start; i < end; ++i) {
			final int itemKey = data[i];
			data[i] = (itemKey ^ ((~itemKey >> Float.SIZE - 1) | (1 << Float.SIZE - 1)));
		}
	}

	public static void encode(long[] data, int start, int end) {
		for (int i = start; i < end; ++i) {
			final long itemKey = data[i];
			data[i] = (itemKey ^ ((itemKey >> Double.SIZE - 1) | (1L << Double.SIZE - 1)));
		}
	}

	public static void decode(long[] data, int start, int end) {
		for (int i = start; i < end; ++i) {
			final long itemKey = data[i];
			data[i] = (itemKey ^ ((~itemKey >> Double.SIZE - 1) | (1L << Double.SIZE - 1)));
		}
	}

	public static int[] radixSortEncodedFloatsWithPermutationTrace(final int[] array, int size) {
		encode(array, 0, size);
		//TODO: move NaNs to end of array and exclude them...
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortUnsignedWithPermutationTrace(array, 0, size, Integer.SIZE - BITS_PER_LEVEL, permutationVector);
		decode(array, 0, size);
		return permutationVector;
	}

	public static int[] radixSortEncodedDoublesWithPermutationTrace(final long[] array, int size) {
		encode(array, 0, size);
		//TODO: move NaNs to end of array and exclude them...
		final int[] permutationVector = generateFreshPermutationVector(size);
		radixSortUnsignedWithPermutationTrace(array, 0, size, Long.SIZE - BITS_PER_LEVEL, permutationVector);
		decode(array, 0, size);
		return permutationVector;
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 * @param permutationVector
	 */
	private static void radixSortUnsignedWithPermutationTrace(final int[] array, final int startOffset, final int end, final int bitShift, final int[] permutationVector) {
		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (array[i] >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortUnsignedWithPermutationTrace(array, startOffset, end, bitShift - BITS_PER_LEVEL, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curPerm = permutationVector[bucketPointers[i]];
				int curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (curItem >> bitShift);
				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpPerm = permutationVector[bucketPointers[itemKey]];
					int tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curPerm;
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					curPerm = tmpPerm;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (curItem >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curPerm;
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortUnsignedWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift, permutationVector);
				} else if (bucketSize > 1) {
					insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector);
				}
			}
		}
	}

	/**
	 * @param array
	 * @param startOffset
	 * @param end
	 * @param bitShift
	 * @param permutationVector
	 */
	private static void radixSortUnsignedWithPermutationTrace(final long[] array, final int startOffset, final int end, final int bitShift, final int[] permutationVector) {
		//this array will initially contain the start offsets for all buckets in this iteration.
		//later on used to remember the current insert position when filling the bucket
		final int[] bucketCounts = new int[NODE_SIZE];
		//this array will contain the end offsets for all buckets in this iteration
		final int[] bucketPointers = new int[NODE_SIZE];

		//determine the number of items in each bucket by counting
		for (int i = end; --i >= startOffset; ) {
			++bucketCounts[MASK & (int) (array[i] >> bitShift)];
		}

		for (int i = 0; i < NODE_SIZE; ++i) {
			int currentCount = bucketCounts[i];
			if (currentCount > 0) {
				if (currentCount == end - startOffset) {
					if (bitShift > 0) {
						radixSortUnsignedWithPermutationTrace(array, startOffset, end, bitShift - BITS_PER_LEVEL, permutationVector);
					}
					return;
				} else {
					break;
				}
			}
		}

		//set the bucketPointers as start indexes and bucketCounts as end indexes for every bucket:
		//first bucket starts on the start offset of the (sub)array
		bucketPointers[0] = startOffset;
		//first bucket ends on start + count.
		bucketCounts[0] += startOffset;
		//now calculate other bucket start and end offsets...
		for (int i = 1; i < NODE_SIZE; ++i) {
			//bucket start = end of predecessor bucket
			bucketPointers[i] = bucketCounts[i - 1];
			//bucket end -> count for bucket + sum of counts of all predecessor buckets
			bucketCounts[i] += bucketCounts[i - 1];
		}

		//for all buckets:
		for (int i = NODE_SIZE; --i >= 0; ) {
			//while current insert position in bucket is not the end...
			while (bucketPointers[i] != bucketCounts[i]) {
				//fetch the current item to distribute into the right bucket
				int curPerm = permutationVector[bucketPointers[i]];
				long curItem = array[bucketPointers[i]];
				//calculate the bucket number as part of some bits of the item
				int itemKey = MASK & (int) (curItem >> bitShift);

				//while we are not in the right bucket for the current item...
				while (i != itemKey) {
					//take the item that currently occupies the slot of our curItem
					int tmpPerm = permutationVector[bucketPointers[itemKey]];
					long tmpItem = array[bucketPointers[itemKey]];
					//insert the curItem in that slot and move the pointer
					permutationVector[bucketPointers[itemKey]] = curPerm;
					array[bucketPointers[itemKey]++] = curItem;
					//make the removed item to curItem
					curItem = tmpItem;
					curPerm = tmpPerm;
					//calculate the key for this item and continue until we find an item belonging to bucket i
					itemKey = MASK & (int) (curItem >> bitShift);
				}
				//i is the right bucket for curItem, so we just insert and move the pointer
				permutationVector[bucketPointers[i]] = curPerm;
				array[bucketPointers[i]++] = curItem;
			}
		}

		//if we are not in the final iteration...
		if (bitShift > 0) {
			//calculate bitshift for the next iteration
			final int nextIterationBitShift = bitShift - BITS_PER_LEVEL;
			//for all buckets...
			for (int i = NODE_SIZE; --i >= 0; ) {
				final int currentPointer = bucketPointers[i];
				//calculate the bucket size
				final int bucketSize = i > 0 ? currentPointer - bucketPointers[i - 1] : currentPointer - startOffset;
				//if bucket size is greater then a threshold recurse, if not switch to insertion sort.
				if (bucketSize > INSERTION_SORT_THRESHOLD) {
					radixSortUnsignedWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, nextIterationBitShift, permutationVector);
				} else if (bucketSize > 1) {
					insertionSortWithPermutationTrace(array, currentPointer - bucketSize, currentPointer, permutationVector);
				}
			}
		}
	}
}

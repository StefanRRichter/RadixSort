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

import sun.misc.Unsafe;

import java.lang.reflect.Field;
import java.nio.ByteOrder;
import java.security.AccessController;
import java.security.PrivilegedAction;

/**
 * Some low level utils.
 */
public final class BitUtils {

	private BitUtils() {
		throw new AssertionError();
	}

	static final Unsafe UNSAFE;
	/**
	 * The offset to the first element in a byte array.
	 */
	static final long BYTE_ARRAY_BASE_OFFSET;
	static final boolean NATIVE_BYTE_ORDER;

	static {
		UNSAFE = (Unsafe) AccessController.doPrivileged(
			new PrivilegedAction<Object>() {

				@Override
				public Object run() {
					try {
						Field f = Unsafe.class.getDeclaredField("theUnsafe");
						f.setAccessible(true);
						return f.get(null);
					} catch (NoSuchFieldException e) {
						throw new Error();
					} catch (IllegalAccessException e) {
						throw new Error();
					}
				}
			});

		BYTE_ARRAY_BASE_OFFSET = UNSAFE != null ? UNSAFE.arrayBaseOffset(byte[].class) : 0;

		// sanity check - this should never fail
		if (UNSAFE != null && UNSAFE.arrayIndexScale(byte[].class) != 1) {
			throw new AssertionError();
		}

		NATIVE_BYTE_ORDER = (ByteOrder.nativeOrder() == ByteOrder.BIG_ENDIAN);
	}

	public static long getLong(byte[] b, int off) {
		if (UNSAFE != null) {
			long ret = UNSAFE.getLong(b, BYTE_ARRAY_BASE_OFFSET + off);
			return NATIVE_BYTE_ORDER ? ret : Long.reverseBytes(ret);
		} else {
			return makeLong(b[off++],
				b[off++],
				b[off++],
				b[off++],
				b[off++],
				b[off++],
				b[off++],
				b[off]);
		}
	}

	static private long makeLong(byte b7, byte b6, byte b5, byte b4,
								 byte b3, byte b2, byte b1, byte b0) {
		return ((((long) b7) << 56)
			| (((long) b6 & 0xff) << 48)
			| (((long) b5 & 0xff) << 40)
			| (((long) b4 & 0xff) << 32)
			| (((long) b3 & 0xff) << 24)
			| (((long) b2 & 0xff) << 16)
			| (((long) b1 & 0xff) << 8)
			| (((long) b0 & 0xff)));
	}

	//--------------------------

	/**
	 * @param v
	 * @return
	 */
	public static boolean containsZeroByte(long v) {
		return (((v) - 0x0101010101010101L) & ~(v) & 0x8080808080808080L) != 0;
	}

	/**
	 * @param data
	 * @param startOffset
	 * @return
	 */
	public static long extractBytesIntoLongOld(final byte[] data, final int startOffset) {
		long key = 0L;
		final int endOffset = Math.min(data.length, startOffset + 8);
		for (int i = startOffset; i < endOffset; ++i) {
			key = key << 8 | (long) data[i];
		}
		key = key << ((8 - endOffset + startOffset) * 8);
		return key;
	}

	/**
	 * @param data
	 * @param startOffset
	 * @return
	 */
	public static long extractBytesIntoLong(final byte[] data, final int startOffset) {
		final int remaining = data.length - startOffset;
		final int endOffset;
		if (remaining >= 8) {
			if (UNSAFE != null) {
				return getLong(data, startOffset);
			}
			endOffset = startOffset + 8;
		} else {
			endOffset = data.length;
		}
		long key = 0L;
		for (int i = startOffset; i < endOffset; ++i) {
			key = key << 8 | (long) data[i];
		}
		key = key << ((8 - endOffset + startOffset) * 8);
		return key;
	}
}
